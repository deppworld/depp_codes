#!/usr/bin/env python3
"""
clonality_check.py  --  Are two tumor samples from the same patient clonally related?

Designed for Illumina TSO500 (DRAGEN) output but works with any VCF + BAM pair.
Runs entirely locally; nothing is uploaded anywhere.

Steps performed
  0. Sample-identity check     : genotype concordance at heterozygous SNPs in both BAMs
                                 (rules out a sample swap before anything else)
  1. Variant harmonisation     : parse T1/T2 VCFs, normalise keys, attach TSO500
                                 CombinedVariantOutput annotations (gene, p., gnomAD) if given
  2. Germline / CHIP filtering : normal VCF, population AF, VAF-shape heuristics, CHIP genes
  3. Force-calling             : every somatic variant from either sample is re-counted
                                 directly from BOTH BAMs by pileup (catches variants the
                                 caller missed in the low-purity sample)
  4. Classification            : shared / T1-private / T2-private, hotspot vs private
  5. Cancer-cell fraction      : purity- and copy-number-corrected CCF per variant
  6. Clonality likelihood      : log10 likelihood ratio (same clone vs independent),
                                 simplified Ostrovnaya-style; report Clonality R for formal stats
  7. 9p21 (CDKN2A/B, MTAP) CN  : per-target coverage ratio in both BAMs, breakpoint
                                 comparison, purity-corrected absolute copy number
  8. Allelic imbalance         : mirrored BAF of germline het SNPs per chromosome arm,
                                 compares LOH landscape between samples
  9. Report                    : clonality_report.txt + TSV tables + PNG plots

Usage (minimal):
  python clonality_check.py \
      --t1-vcf T1.vcf.gz --t2-vcf T2.vcf.gz \
      --t1-bam T1.bam    --t2-bam T2.bam \
      --purity1 0.75 --purity2 0.30 \
      --targets TSO500_targets.bed --genome hg19 \
      --out clonality_out

Optional but recommended:
  --normal-vcf blood.vcf.gz            matched germline
  --t1-cvo T1_CombinedVariantOutput.tsv --t2-cvo T2_CombinedVariantOutput.tsv
  --cn T1_T2_gene_cn.tsv               (gene<TAB>cn_T1<TAB>cn_T2) for CCF correction
  --hotspots hotspots.tsv               (gene<TAB>protein_regex<TAB>recurrence_freq)
"""

import argparse
import gzip
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

CHIP_GENES = {"DNMT3A", "TET2", "ASXL1", "JAK2", "SF3B1", "SRSF2", "TP53", "PPM1D",
              "GNB1", "CBL", "U2AF1", "IDH2", "KRAS", "GNAS", "BCOR", "STAG2"}
# NB: TP53/KRAS/IDH2 are listed because low-VAF (<3%) calls with identical VAF in
# both samples can be CHIP; they are only flagged, never auto-removed.

# Built-in hotspot table (gene, regex on p. notation, approximate recurrence
# frequency across all cancers).  Used only if --hotspots not supplied.
DEFAULT_HOTSPOTS = [
    ("TP53",   r"p\.\(?(R175|Y220|G245|R248|R249|R273|R282)", 0.03),
    ("KRAS",   r"p\.\(?(G12|G13|Q61)", 0.10),
    ("NRAS",   r"p\.\(?(G12|G13|Q61)", 0.02),
    ("HRAS",   r"p\.\(?(G12|G13|Q61)", 0.01),
    ("BRAF",   r"p\.\(?V600", 0.05),
    ("PIK3CA", r"p\.\(?(E542|E545|H1047)", 0.05),
    ("CTNNB1", r"p\.\(?(D32|S33|G34|S37|T41|S45)", 0.02),
    ("IDH1",   r"p\.\(?R132", 0.02),
    ("IDH2",   r"p\.\(?(R140|R172)", 0.01),
    ("EGFR",   r"p\.\(?(L858|T790|G719|L861)", 0.02),
    ("AKT1",   r"p\.\(?E17K", 0.01),
    ("GNAS",   r"p\.\(?R201", 0.01),
    ("FGFR3",  r"p\.\(?(R248|S249|G370|Y373)", 0.01),
    ("TERT",   r"promoter|c\.-(124|146)", 0.08),
]

# 9p21 region (hg19 and hg38), 1-based inclusive
REGIONS_9P21 = {
    "hg19": {"window": ("9", 21_700_000, 22_300_000),
             "MTAP":   (21_802_636, 21_865_970),
             "CDKN2A": (21_967_752, 21_995_301),
             "CDKN2B": (22_002_903, 22_009_313)},
    "hg38": {"window": ("9", 21_700_000, 22_300_000),
             "MTAP":   (21_802_636, 21_865_970),
             "CDKN2A": (21_967_752, 21_995_301),
             "CDKN2B": (22_002_903, 22_009_313)},
}
# (CDKN2A/B/MTAP coordinates are identical between hg19 and hg38 on chr9 for
#  practical purposes; the window is chosen wide enough to cover both.)

# Approximate centromere positions (Mb) for arm assignment
CENTROMERE_MB = {"1": 125, "2": 93, "3": 91, "4": 50, "5": 48, "6": 61, "7": 60,
                 "8": 45, "9": 49, "10": 40, "11": 53, "12": 35, "13": 17, "14": 17,
                 "15": 19, "16": 36, "17": 25, "18": 17, "19": 26, "20": 27,
                 "21": 13, "22": 15, "X": 60}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def norm_chrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def bam_chrom(bam, chrom):
    """Return the contig name as spelled in this BAM ('9' or 'chr9')."""
    refs = set(bam.references)
    if chrom in refs:
        return chrom
    if "chr" + chrom in refs:
        return "chr" + chrom
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]
    return None


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def log(msg):
    print(f"[clonality] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# Step 1: parse VCFs
# --------------------------------------------------------------------------- #

def parse_vcf(path, sample_label, min_dp=20, min_alt=3, min_vaf=0.02, pass_only=True):
    """Return DataFrame of variants with key, vaf, dp, alt reads, filter."""
    vf = pysam.VariantFile(path)
    rows = []
    for rec in vf:
        if rec.alts is None:
            continue
        # gVCF non-variant blocks
        if rec.alts == ("<NON_REF>",) or rec.alts[0] is None:
            continue
        filt = ",".join(rec.filter.keys()) if rec.filter.keys() else "PASS"
        if pass_only and filt not in ("PASS", "."):
            continue
        for ai, alt in enumerate(rec.alts):
            if alt in (None, "<NON_REF>", "*"):
                continue
            dp, ad_ref, ad_alt, vaf = None, None, None, None
            if len(vf.header.samples) > 0:
                s = rec.samples[0]
                dp = s.get("DP")
                ad = s.get("AD")
                if ad is not None and len(ad) > ai + 1:
                    ad_ref, ad_alt = ad[0], ad[ai + 1]
                af = s.get("AF")
                if af is not None:
                    vaf = af[ai] if isinstance(af, tuple) else af
            if vaf is None and ad_alt is not None and dp:
                vaf = ad_alt / dp
            if vaf is None and ad_alt is not None and ad_ref is not None and (ad_ref + ad_alt) > 0:
                vaf = ad_alt / (ad_ref + ad_alt)
            if dp is None:
                dp = rec.info.get("DP")
            if vaf is None or dp is None:
                continue
            if dp < min_dp or vaf < min_vaf:
                continue
            if ad_alt is not None and ad_alt < min_alt:
                continue
            rows.append({
                "key": f"{norm_chrom(rec.chrom)}:{rec.pos}:{rec.ref}:{alt}",
                "chrom": norm_chrom(rec.chrom), "pos": rec.pos, "ref": rec.ref, "alt": alt,
                f"vaf_{sample_label}_called": float(vaf),
                f"dp_{sample_label}_called": int(dp),
                f"filter_{sample_label}": filt,
            })
    df = pd.DataFrame(rows)
    log(f"{sample_label}: {len(df)} variants parsed from {os.path.basename(path)}")
    return df


def parse_cvo(path):
    """Parse the [Small Variants] block of a TSO500 CombinedVariantOutput.tsv."""
    if path is None:
        return pd.DataFrame(columns=["key", "gene", "pdot", "cdot", "consequence"])
    rows, in_block, header = [], False, None
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("[Small Variants]"):
                in_block, header = True, None
                continue
            if in_block and line.startswith("["):
                break
            if not in_block or not line.strip():
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            d = dict(zip(header, parts))
            try:
                key = (f"{norm_chrom(d['Chromosome'])}:{int(d['Genomic Position'])}:"
                       f"{d['Reference Call']}:{d['Alternative Call']}")
            except (KeyError, ValueError):
                continue
            rows.append({"key": key, "gene": d.get("Gene", ""),
                         "pdot": d.get("P-Dot Notation", ""),
                         "cdot": d.get("C-Dot Notation", ""),
                         "consequence": d.get("Consequence(s)", "")})
    df = pd.DataFrame(rows).drop_duplicates("key")
    log(f"CVO annotations: {len(df)} small variants from {os.path.basename(path)}")
    return df


def parse_tmb_trace(path, min_pop_alleles=5):
    """Parse a TSO500 TMB trace TSV (<sample>_TMB_Trace.tsv or <sample>.dna.tmb.trace.tsv).

    Returns (germline_keys, somatic_keys).  A variant is germline if any
    'GermlineFilter*' column is True, or if population allele counts
    (gnomAD exome/genome, 1000G) reach min_pop_alleles.
    """
    germ, som = set(), set()
    if path is None:
        return germ, som
    with open_text(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                hl = [h.lower() for h in header]
                gcols = [i for i, h in enumerate(hl) if "germlinefilter" in h]
                pcols = [i for i, h in enumerate(hl) if "allelecount" in h and ("gnomad" in h or "1000" in h)]
                try:
                    ic = next(i for i, h in enumerate(hl) if h in ("chromosome", "chrom", "chr"))
                    ip = next(i for i, h in enumerate(hl) if h in ("position", "pos", "genomic position"))
                    ir = next(i for i, h in enumerate(hl) if h in ("refcall", "ref", "reference call"))
                    ia = next(i for i, h in enumerate(hl) if h in ("altcall", "alt", "alternative call"))
                except StopIteration:
                    log(f"TMB trace {path}: could not find chrom/pos/ref/alt columns; ignored")
                    return germ, som
                continue
            if len(parts) < len(header):
                continue
            key = f"{norm_chrom(parts[ic])}:{parts[ip]}:{parts[ir]}:{parts[ia]}"
            is_germ = any(parts[i].strip().lower() in ("true", "1", "yes") for i in gcols)
            for i in pcols:
                try:
                    if float(parts[i]) >= min_pop_alleles:
                        is_germ = True
                except ValueError:
                    pass
            (germ if is_germ else som).add(key)
    log(f"TMB trace {os.path.basename(path)}: {len(germ)} germline-flagged, {len(som)} somatic-classified")
    return germ, som


def load_hotspots(path):
    if path is None:
        return DEFAULT_HOTSPOTS
    hs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            g, rx, f = line.rstrip("\n").split("\t")[:3]
            hs.append((g, rx, float(f)))
    return hs


def hotspot_freq(gene, pdot, cdot, hotspots):
    for g, rx, f in hotspots:
        if gene == g and (re.search(rx, str(pdot) or "") or re.search(rx, str(cdot) or "")):
            return f
    return None


# --------------------------------------------------------------------------- #
# Step 3: force-calling from BAM
# --------------------------------------------------------------------------- #

def count_alleles(bam, chrom, pos, ref, alt, min_bq=20, min_mq=20, max_depth=100000):
    """Count ref/alt supporting reads at a variant by pileup.

    SNV  : base match at pos
    DEL  : ref longer than alt; count reads whose pileup 'indel' == -(len diff) at anchor
    INS  : alt longer than ref; count reads with indel == +len diff at anchor
    Returns (ref_count, alt_count, depth).
    """
    c = bam_chrom(bam, chrom)
    if c is None:
        return 0, 0, 0
    is_snv = len(ref) == 1 and len(alt) == 1
    is_del = len(ref) > len(alt)
    is_ins = len(alt) > len(ref)
    indel_len = len(alt) - len(ref)
    n_ref = n_alt = depth = 0
    for col in bam.pileup(c, pos - 1, pos, truncate=True, min_base_quality=min_bq,
                          min_mapping_quality=min_mq, max_depth=max_depth,
                          ignore_overlaps=True, stepper="samtools"):
        if col.reference_pos != pos - 1:
            continue
        for pr in col.pileups:
            depth += 1
            if pr.is_refskip:
                continue
            if is_snv:
                if pr.is_del or pr.query_position is None:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                if b == alt.upper():
                    n_alt += 1
                elif b == ref.upper():
                    n_ref += 1
            elif is_del or is_ins:
                if pr.indel == indel_len and indel_len != 0:
                    n_alt += 1
                elif pr.indel == 0 and not pr.is_del:
                    n_ref += 1
            else:  # MNV / complex: approximate by first base
                if pr.is_del or pr.query_position is None:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                if b == alt[0].upper() and alt[0] != ref[0]:
                    n_alt += 1
                elif b == ref[0].upper():
                    n_ref += 1
    return n_ref, n_alt, depth


def force_call(df, bam1, bam2, min_bq, min_mq):
    out = defaultdict(list)
    for _, r in df.iterrows():
        for lab, bam in (("T1", bam1), ("T2", bam2)):
            nr, na, dp = count_alleles(bam, r.chrom, r.pos, r.ref, r.alt, min_bq, min_mq)
            out[f"ref_{lab}"].append(nr)
            out[f"alt_{lab}"].append(na)
            out[f"dp_{lab}"].append(dp)
            out[f"vaf_{lab}"].append(na / (nr + na) if (nr + na) > 0 else 0.0)
    for k, v in out.items():
        df[k] = v
    return df


# --------------------------------------------------------------------------- #
# Step 5: cancer-cell fraction
# --------------------------------------------------------------------------- #

def ccf(vaf, purity, cn_tumor=2, multiplicity=1):
    """CCF = VAF * (purity*CN_t + 2*(1-purity)) / (purity * multiplicity). Clipped to [0, 1.5]."""
    if purity <= 0 or vaf is None:
        return float("nan")
    val = vaf * (purity * cn_tumor + 2 * (1 - purity)) / (purity * multiplicity)
    return float(min(val, 1.5))


def load_cn(path):
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = [c.lower() for c in df.columns]
    return {r["gene"]: (float(r["cn_t1"]), float(r["cn_t2"])) for _, r in df.iterrows()}


# --------------------------------------------------------------------------- #
# Step 6: clonality likelihood ratio
# --------------------------------------------------------------------------- #

def clonality_lr(shared, panel_size_bp, n_som_t1, n_som_t2):
    """
    Simplified likelihood ratio, in the spirit of Ostrovnaya et al. (Clonality R pkg).

    For a shared PRIVATE (non-hotspot) variant, P(chance match | independent) is
    approximated by the probability that the second tumor independently mutates the
    same base to the same allele:  p ~ (n_somatic_in_other_tumor / panel_bp) / 3.
    For a shared HOTSPOT the probability is the hotspot recurrence frequency.
    LR contribution per shared variant = 1 / p.  Variants private to one sample do not
    penalise (metastases and subclones lose/gain variants), so this is a
    one-sided, conservative-in-the-other-direction statistic: interpret log10 LR >= 3
    as strong support for a common clonal origin, 1-3 as moderate, < 1 as
    uninformative.  Use Clonality R (clonalityAnalysis / LRtesting) for a formal
    test with a proper reference set.
    """
    per_site = max(n_som_t1, n_som_t2, 1) / max(panel_size_bp, 1) / 3.0
    shared = shared.copy()
    shared["event_id"] = collapse_events(shared)
    log10lr = 0.0
    contribs = []
    seen = set()
    for _, r in shared.iterrows():
        p = r["hotspot_freq"] if not pd.isna(r["hotspot_freq"]) else per_site
        p = min(max(p, 1e-9), 0.5)
        c = -math.log10(p)
        if r["event_id"] in seen:      # same complex indel written as several VCF records
            c = 0.0
        seen.add(r["event_id"])
        contribs.append(c)
        log10lr += c
    shared["log10_LR_contribution"] = contribs
    return log10lr, shared


def collapse_events(df, window=30):
    """Assign one event id to VCF records that are within `window` bp on the same
    chromosome and involve an indel or MNV (DRAGEN writes complex indels as several
    overlapping records).  SNVs are never merged with each other."""
    ids = {}
    n = 0
    rows = sorted(df.itertuples(), key=lambda r: (r.key.split(":")[0], int(r.key.split(":")[1])))
    last_chrom, last_pos, last_id, last_indel = None, -10**9, None, False
    for r in rows:
        chrom, pos, ref, alt = r.key.split(":")[:4]
        pos = int(pos)
        indel = len(ref) != 1 or len(alt) != 1
        if chrom == last_chrom and pos - last_pos <= window and indel and last_indel:
            ids[r.key] = last_id
        else:
            n += 1
            ids[r.key] = f"E{n}"
            last_id = f"E{n}"
        last_chrom, last_pos, last_indel = chrom, pos, indel
    return [ids[k] for k in df.key]


# --------------------------------------------------------------------------- #
# Step 7: 9p21 coverage
# --------------------------------------------------------------------------- #

def read_bed(path):
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            rows.append({"chrom": norm_chrom(p[0]), "start": int(p[1]), "end": int(p[2]),
                         "name": p[3] if len(p) > 3 else f"{p[0]}:{p[1]}-{p[2]}"})
    return pd.DataFrame(rows)


def target_depths(bam, targets, min_mq=20):
    vals = []
    for _, t in targets.iterrows():
        c = bam_chrom(bam, t.chrom)
        if c is None:
            vals.append(np.nan)
            continue
        cov = bam.count_coverage(c, t.start, t.end, quality_threshold=0,
                                 read_callback=lambda r: r.mapping_quality >= min_mq and not r.is_duplicate)
        depth = np.sum(cov, axis=0)
        vals.append(float(np.median(depth)) if len(depth) else np.nan)
    return np.array(vals)


def cn_from_ratio(ratio, purity):
    """Absolute tumor copy number implied by a coverage ratio at given purity
    (ratio normalised to diploid autosomes)."""
    if purity <= 0:
        return float("nan")
    return (2 * ratio - 2 * (1 - purity)) / purity


def analyse_9p21(bam1, bam2, targets, genome, purity1, purity2, outdir, min_mq):
    reg = REGIONS_9P21[genome]
    wchrom, wstart, wend = reg["window"]
    win = targets[(targets.chrom == wchrom) & (targets.end > wstart) & (targets.start < wend)].copy()
    if win.empty:
        log("No target regions found in the 9p21 window; skipping CN analysis "
            "(check --targets and --genome).")
        return None
    # normalisation set: all autosomal targets outside chr9
    norm = targets[(targets.chrom != "9") & (targets.chrom != "X") & (targets.chrom != "Y")]
    if len(norm) > 400:  # subsample for speed
        norm = norm.sample(400, random_state=1)
    log(f"9p21: {len(win)} targets in window, {len(norm)} normalisation targets")
    res = win.reset_index(drop=True)
    for lab, bam in (("T1", bam1), ("T2", bam2)):
        d_win = target_depths(bam, win, min_mq)
        d_norm = target_depths(bam, norm, min_mq)
        med = np.nanmedian(d_norm)
        res[f"depth_{lab}"] = d_win
        res[f"ratio_{lab}"] = d_win / med if med and med > 0 else np.nan
    res["gene"] = ""
    for g in ("MTAP", "CDKN2A", "CDKN2B"):
        s, e = reg[g]
        res.loc[(res.end > s) & (res.start < e), "gene"] = g
    res["mid"] = (res.start + res.end) / 2
    res = res.sort_values("mid")

    # purity-corrected absolute CN over CDKN2A targets
    summary = {}
    for lab, pur in (("T1", purity1), ("T2", purity2)):
        for g in ("MTAP", "CDKN2A", "CDKN2B"):
            r = res.loc[res.gene == g, f"ratio_{lab}"].median()
            summary[(lab, g)] = (r, cn_from_ratio(r, pur) if not np.isnan(r) else np.nan)
    # crude breakpoint estimate: contiguous run of targets with ratio < threshold
    bp = {}
    for lab, pur in (("T1", purity1), ("T2", purity2)):
        # expected ratio for CN=0 at this purity
        thr = (1 - pur) + 0.15
        low = res[res[f"ratio_{lab}"] < thr]
        bp[lab] = (int(low.start.min()), int(low.end.max())) if not low.empty else None

    res.to_csv(os.path.join(outdir, "cn_9p21_targets.tsv"), sep="\t", index=False)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4))
        for lab, col in (("T1", "#1f77b4"), ("T2", "#d62728")):
            ax.plot(res.mid / 1e6, res[f"ratio_{lab}"], "o-", ms=4, label=lab, color=col)
        for g in ("MTAP", "CDKN2A", "CDKN2B"):
            s, e = reg[g]
            ax.axvspan(s / 1e6, e / 1e6, alpha=0.15, color="grey")
            ax.text((s + e) / 2e6, 1.45, g, ha="center", fontsize=8)
        ax.axhline(1.0, ls="--", color="k", lw=0.8)
        ax.set_ylim(0, 1.6)
        ax.set_xlabel(f"chr9 position (Mb, {genome})")
        ax.set_ylabel("coverage ratio (vs off-chr9 targets)")
        ax.set_title("9p21 per-target coverage: T1 vs T2")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "cn_9p21_profile.png"), dpi=150)
        plt.close(fig)
    except Exception as e:  # plotting must never kill the run
        log(f"plot skipped: {e}")
    return res, summary, bp


# --------------------------------------------------------------------------- #
# Step 0 / 8: sample identity and allelic imbalance from het SNPs
# --------------------------------------------------------------------------- #

def het_snp_candidates(vcf_paths, min_dp=50):
    """Collect SNVs that look germline heterozygous (VAF 0.35-0.65) in any VCF."""
    cand = {}
    for p in vcf_paths:
        if p is None:
            continue
        vf = pysam.VariantFile(p)
        for rec in vf:
            if rec.alts is None or len(rec.ref) != 1:
                continue
            for ai, alt in enumerate(rec.alts):
                if alt is None or len(alt) != 1 or alt == "<NON_REF>":
                    continue
                s = rec.samples[0] if len(vf.header.samples) else None
                if s is None:
                    continue
                dp = s.get("DP")
                af = s.get("AF")
                if af is None:
                    ad = s.get("AD")
                    if ad and dp:
                        af = ad[ai + 1] / dp
                    else:
                        continue
                af = af[ai] if isinstance(af, tuple) else af
                if dp and dp >= min_dp and 0.35 <= af <= 0.65:
                    cand[f"{norm_chrom(rec.chrom)}:{rec.pos}:{rec.ref}:{alt}"] = \
                        (norm_chrom(rec.chrom), rec.pos, rec.ref, alt)
    return cand


def identity_and_baf(cand, bam1, bam2, outdir, min_bq, min_mq, max_snps=3000):
    keys = list(cand.keys())[:max_snps]
    rows = []
    for k in keys:
        chrom, pos, ref, alt = cand[k]
        r1 = count_alleles(bam1, chrom, pos, ref, alt, min_bq, min_mq)
        r2 = count_alleles(bam2, chrom, pos, ref, alt, min_bq, min_mq)
        rows.append({"key": k, "chrom": chrom, "pos": pos,
                     "vaf_T1": r1[1] / (r1[0] + r1[1]) if (r1[0] + r1[1]) >= 20 else np.nan,
                     "vaf_T2": r2[1] / (r2[0] + r2[1]) if (r2[0] + r2[1]) >= 20 else np.nan})
    df = pd.DataFrame(rows).dropna()
    if df.empty:
        return None

    def gt(v):
        return 0 if v < 0.15 else (2 if v > 0.85 else 1)

    df["gt_T1"] = df.vaf_T1.map(gt)
    df["gt_T2"] = df.vaf_T2.map(gt)
    # identity: a SNP het in one sample must be at least present (>0.15) in the other.
    # Complete absence (0/0 vs 1/1 or 1/1 vs 0/0) is the sample-swap signature.
    discordant = ((df.gt_T1 == 0) & (df.gt_T2 == 2)) | ((df.gt_T1 == 2) & (df.gt_T2 == 0))
    n_hom_disc = int(discordant.sum())
    absent_disc = ((df.gt_T1 == 0) & (df.gt_T2 > 0)) | ((df.gt_T2 == 0) & (df.gt_T1 > 0))
    frac_absent = float(absent_disc.mean())
    # allelic imbalance per arm (mirrored BAF)
    df["arm"] = [f"{c}{'p' if p / 1e6 < CENTROMERE_MB.get(c, 1e9) else 'q'}"
                 for c, p in zip(df.chrom, df.pos)]
    df["mbaf_T1"] = (df.vaf_T1 - 0.5).abs()
    df["mbaf_T2"] = (df.vaf_T2 - 0.5).abs()
    arm = df.groupby("arm").agg(n=("key", "size"), mbaf_T1=("mbaf_T1", "median"),
                                mbaf_T2=("mbaf_T2", "median")).reset_index()
    arm = arm[arm.n >= 5]
    df.to_csv(os.path.join(outdir, "het_snps_T1_T2.tsv"), sep="\t", index=False)
    arm.to_csv(os.path.join(outdir, "allelic_imbalance_by_arm.tsv"), sep="\t", index=False)
    corr = float(np.corrcoef(arm.mbaf_T1, arm.mbaf_T2)[0, 1]) if len(arm) >= 4 else float("nan")
    return {"n_snps": len(df), "n_hom_discordant": n_hom_disc, "frac_absent_discordant": frac_absent,
            "arm_table": arm, "arm_corr": corr}


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1-vcf", required=True)
    ap.add_argument("--t2-vcf", required=True)
    ap.add_argument("--t1-bam", required=True)
    ap.add_argument("--t2-bam", required=True)
    ap.add_argument("--normal-vcf", default=None, help="matched germline VCF (blood/normal)")
    ap.add_argument("--t1-cvo", default=None, help="TSO500 CombinedVariantOutput.tsv for T1")
    ap.add_argument("--t2-cvo", default=None, help="TSO500 CombinedVariantOutput.tsv for T2")
    ap.add_argument("--t1-tmb-trace", default=None, help="TSO500 TMB trace TSV for T1 (germline flags)")
    ap.add_argument("--t2-tmb-trace", default=None, help="TSO500 TMB trace TSV for T2 (germline flags)")
    ap.add_argument("--germline-purity-margin", type=float, default=0.10,
                    help="VAF above (purity + margin) in either sample => germline (set 1.0 to disable)")
    ap.add_argument("--purity1", type=float, required=True, help="tumor purity T1 (0-1), from pathology or ploidy tool")
    ap.add_argument("--purity2", type=float, required=True, help="tumor purity T2 (0-1)")
    ap.add_argument("--targets", default=None, help="panel target BED (needed for 9p21 CN step)")
    ap.add_argument("--genome", default="hg19", choices=["hg19", "hg38"])
    ap.add_argument("--cn", default=None, help="TSV gene<TAB>cn_T1<TAB>cn_T2 for CCF correction")
    ap.add_argument("--hotspots", default=None, help="TSV gene<TAB>protein_regex<TAB>recurrence_freq")
    ap.add_argument("--panel-size-bp", type=float, default=1.94e6, help="callable panel size (TSO500 ~1.94 Mb)")
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--min-alt", type=int, default=3)
    ap.add_argument("--min-vaf", type=float, default=0.02, help="VAF floor for a caller call to be considered")
    ap.add_argument("--force-min-alt", type=int, default=3, help="alt reads needed to call a variant present by force-calling")
    ap.add_argument("--force-min-vaf", type=float, default=0.005)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--germline-vaf-band", type=float, nargs=2, default=(0.35, 0.65),
                    help="VAF band in BOTH samples that flags a variant as likely germline when no normal is given")
    ap.add_argument("--no-pass-only", action="store_true", help="also consider non-PASS calls")
    ap.add_argument("--out", default="clonality_out")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    bam1 = pysam.AlignmentFile(a.t1_bam, "rb")
    bam2 = pysam.AlignmentFile(a.t2_bam, "rb")
    hotspots = load_hotspots(a.hotspots)
    cn_map = load_cn(a.cn)
    report = []

    def R(s=""):
        report.append(s)
        print(s)

    R("=" * 78)
    R("CLONALITY ASSESSMENT  T1 vs T2")
    R("=" * 78)
    R(f"T1: {a.t1_bam}  purity={a.purity1}")
    R(f"T2: {a.t2_bam}  purity={a.purity2}")
    R(f"Normal VCF: {a.normal_vcf or 'NONE (germline inferred heuristically - weaker)'}")
    R()

    # ---- Step 0: sample identity + allelic imbalance -------------------------
    log("Step 0/8: sample identity check from heterozygous SNPs")
    cand = het_snp_candidates([a.t1_vcf, a.t2_vcf, a.normal_vcf])
    ident = identity_and_baf(cand, bam1, bam2, a.out, a.min_bq, a.min_mq) if cand else None
    R("-" * 78)
    R("STEP 0  SAMPLE IDENTITY (same patient?)")
    if ident is None:
        R("  Not enough heterozygous SNPs found in VCFs (TSO500 gVCF or a --normal-vcf helps).")
    else:
        R(f"  het SNPs evaluated in both BAMs : {ident['n_snps']}")
        R(f"  hom-ref vs hom-alt discordances : {ident['n_hom_discordant']}")
        R(f"  fraction absent in one sample   : {ident['frac_absent_discordant']:.3f}")
        if ident["n_hom_discordant"] > max(2, 0.02 * ident["n_snps"]) or ident["frac_absent_discordant"] > 0.20:
            R("  >>> WARNING: genotype discordance is high. SAMPLE SWAP / DIFFERENT PATIENT "
              "must be excluded before any clonality interpretation. <<<")
        else:
            R("  Genotypes concordant: consistent with the same patient.")
    R()

    # ---- Step 1: parse -------------------------------------------------------
    log("Step 1/8: parsing VCFs")
    t1 = parse_vcf(a.t1_vcf, "T1", a.min_dp, a.min_alt, a.min_vaf, not a.no_pass_only)
    t2 = parse_vcf(a.t2_vcf, "T2", a.min_dp, a.min_alt, a.min_vaf, not a.no_pass_only)
    if t1.empty and t2.empty:
        R("No variants passed filters in either VCF. Nothing to compare.")
        sys.exit(1)
    allv = pd.concat([t1, t2]).drop_duplicates("key").reset_index(drop=True)
    for lab, df in (("T1", t1), ("T2", t2)):
        allv = allv.merge(df[["key", f"vaf_{lab}_called", f"dp_{lab}_called", f"filter_{lab}"]],
                          on="key", how="left", suffixes=("", "_dup"))
        for c in list(allv.columns):
            if c.endswith("_dup"):
                base = c[:-4]
                allv[base] = allv[base].fillna(allv[c])
                allv.drop(columns=c, inplace=True)
    cvo = pd.concat([parse_cvo(a.t1_cvo), parse_cvo(a.t2_cvo)]).drop_duplicates("key")
    allv = allv.merge(cvo, on="key", how="left")
    for c in ("gene", "pdot", "cdot", "consequence"):
        allv[c] = allv[c].fillna("")

    # ---- Step 2: germline / CHIP flags ---------------------------------------
    log("Step 2/8: germline / CHIP filtering")
    normal_keys = set()
    if a.normal_vcf:
        n = parse_vcf(a.normal_vcf, "N", min_dp=10, min_alt=2, min_vaf=0.15, pass_only=False)
        normal_keys = set(n.key)
    allv["in_normal"] = allv.key.isin(normal_keys)
    g1, s1 = parse_tmb_trace(a.t1_tmb_trace)
    g2, s2 = parse_tmb_trace(a.t2_tmb_trace)
    allv["tso500_germline"] = allv.key.isin(g1 | g2)

    # ---- Step 3: force-call ----------------------------------------------------
    log(f"Step 3/8: force-calling {len(allv)} variants in both BAMs")
    allv = force_call(allv, bam1, bam2, a.min_bq, a.min_mq)
    lo, hi = a.germline_vaf_band
    # Purity-aware rule: a somatic VAF cannot exceed the sample's purity (except with
    # mutant-allele amplification), whereas germline VAFs are purity-independent.
    # A variant above purity+margin in EITHER sample is treated as germline unless it
    # is a known hotspot (hotspot check happens below, so we exclude those afterwards).
    m = a.germline_purity_margin
    over_purity = (allv.vaf_T1 > a.purity1 + m) | (allv.vaf_T2 > a.purity2 + m)
    allv["germline_like"] = allv.in_normal | allv.tso500_germline | (
        allv.vaf_T1.between(lo, hi) & allv.vaf_T2.between(lo, hi)) | (
        (allv.vaf_T1 > 0.85) & (allv.vaf_T2 > 0.85)) | over_purity
    allv["hotspot_freq"] = [hotspot_freq(g, p, c, hotspots) for g, p, c in zip(allv.gene, allv.pdot, allv.cdot)]
    allv["hotspot_freq"] = allv.hotspot_freq.astype(float)
    is_hot = allv.hotspot_freq.notna()
    # hotspot drivers with LOH / amplification may legitimately exceed purity: keep them
    allv.loc[is_hot & over_purity & ~allv.in_normal & ~allv.tso500_germline, "germline_like"] = False
    allv["chip_flag"] = allv.gene.isin(CHIP_GENES) & (allv.vaf_T1 < 0.05) & (allv.vaf_T2 < 0.05) & \
        ((allv.vaf_T1 - allv.vaf_T2).abs() < 0.02) & (allv.alt_T1 >= 3) & (allv.alt_T2 >= 3)

    # ---- Step 4: classification ---------------------------------------------
    log("Step 4/8: classification")
    present1 = (allv.alt_T1 >= a.force_min_alt) & (allv.vaf_T1 >= a.force_min_vaf)
    present2 = (allv.alt_T2 >= a.force_min_alt) & (allv.vaf_T2 >= a.force_min_vaf)
    allv["status"] = np.select(
        [allv.germline_like, present1 & present2, present1 & ~present2, ~present1 & present2],
        ["GERMLINE_LIKE", "SHARED", "T1_PRIVATE", "T2_PRIVATE"], default="UNSUPPORTED")
    allv["variant_class"] = np.where(allv.hotspot_freq.notna(), "HOTSPOT", "PRIVATE")
    # 'present in T2 only by force-calling' = caller missed it (typical for low purity)
    allv["rescued_in_T2"] = allv["vaf_T2_called"].isna() & present2 & present1
    allv["rescued_in_T1"] = allv["vaf_T1_called"].isna() & present1 & present2

    # ---- Step 5: CCF -----------------------------------------------------------
    log("Step 5/8: cancer-cell fraction")
    def cn_for(gene, idx):
        return cn_map.get(gene, (2.0, 2.0))[idx] if gene else 2.0
    allv["ccf_T1"] = [ccf(v, a.purity1, cn_for(g, 0)) for v, g in zip(allv.vaf_T1, allv.gene)]
    allv["ccf_T2"] = [ccf(v, a.purity2, cn_for(g, 1)) for v, g in zip(allv.vaf_T2, allv.gene)]
    allv["clonal_T1"] = allv.ccf_T1 >= 0.7
    allv["clonal_T2"] = allv.ccf_T2 >= 0.7

    cols = ["key", "gene", "pdot", "consequence", "status", "variant_class", "hotspot_freq",
            "vaf_T1_called", "vaf_T2_called", "dp_T1", "alt_T1", "vaf_T1", "dp_T2", "alt_T2", "vaf_T2",
            "ccf_T1", "ccf_T2", "clonal_T1", "clonal_T2", "rescued_in_T1", "rescued_in_T2",
            "in_normal", "tso500_germline", "germline_like", "chip_flag", "filter_T1", "filter_T2"]
    allv = allv[cols].sort_values(["status", "gene", "key"])
    allv.to_csv(os.path.join(a.out, "variants_T1_T2.tsv"), sep="\t", index=False, float_format="%.4f")

    somatic = allv[allv.status.isin(["SHARED", "T1_PRIVATE", "T2_PRIVATE"]) & ~allv.chip_flag]
    shared = somatic[somatic.status == "SHARED"]
    p1 = somatic[somatic.status == "T1_PRIVATE"]
    p2 = somatic[somatic.status == "T2_PRIVATE"]

    R("-" * 78)
    R("STEP 1-4  VARIANT OVERLAP (after force-calling both BAMs)")
    R(f"  germline-like (excluded)  : {int((allv.status == 'GERMLINE_LIKE').sum())}  "
      f"(TSO500-flagged {int(allv.tso500_germline.sum())}, in normal {int(allv.in_normal.sum())}, "
      f"VAF>purity+{a.germline_purity_margin:.2f} or germline VAF band: rest)")
    R(f"  possible CHIP (excluded)  : {int(allv.chip_flag.sum())}")
    R(f"  SHARED somatic            : {len(shared)}  "
      f"(private/non-hotspot: {int((shared.variant_class == 'PRIVATE').sum())}, "
      f"hotspot: {int((shared.variant_class == 'HOTSPOT').sum())})")
    R(f"     of which rescued in T2 by force-calling (missed by caller): {int(shared.rescued_in_T2.sum())}")
    R(f"     of which rescued in T1 by force-calling                  : {int(shared.rescued_in_T1.sum())}")
    R(f"  T1-private somatic        : {len(p1)}")
    R(f"  T2-private somatic        : {len(p2)}")
    R()
    if not shared.empty:
        R("  Shared variants:")
        R("  {:<28s} {:<9s} {:<18s} {:>7s} {:>7s} {:>6s} {:>6s}  {}".format(
            "key", "gene", "p.", "VAF_T1", "VAF_T2", "CCF_T1", "CCF_T2", "class"))
        for _, r in shared.iterrows():
            R("  {:<28s} {:<9s} {:<18s} {:>7.3f} {:>7.3f} {:>6.2f} {:>6.2f}  {}{}".format(
                r.key[:28], r.gene[:9], str(r.pdot)[:18], r.vaf_T1, r.vaf_T2, r.ccf_T1, r.ccf_T2,
                r.variant_class, " (rescued)" if r.rescued_in_T2 or r.rescued_in_T1 else ""))
        R()
    for lab, dfp in (("T1", p1), ("T2", p2)):
        if not dfp.empty:
            R(f"  {lab}-private variants:")
            for _, r in dfp.iterrows():
                R("  {:<28s} {:<9s} {:<18s} VAF_T1={:.3f} VAF_T2={:.3f} CCF_{}={:.2f}".format(
                    r.key[:28], r.gene[:9], str(r.pdot)[:18], r.vaf_T1, r.vaf_T2, lab,
                    r.ccf_T1 if lab == "T1" else r.ccf_T2))
            R()

    # ---- Step 6: LR ------------------------------------------------------------
    log("Step 6/8: clonality likelihood ratio")
    n1 = len(shared) + len(p1)
    n2 = len(shared) + len(p2)
    log10lr, shared_lr = clonality_lr(shared, a.panel_size_bp, n1, n2)
    shared_lr.to_csv(os.path.join(a.out, "shared_variants_LR.tsv"), sep="\t", index=False, float_format="%.4f")
    n_events = shared_lr.event_id.nunique() if not shared_lr.empty else 0
    n_priv_shared = int(shared_lr.drop_duplicates("event_id").variant_class.eq("PRIVATE").sum()) if n_events else 0
    R("-" * 78)
    R("STEP 6  CLONALITY LIKELIHOOD (same clone vs independent tumors)")
    R(f"  distinct shared events (overlapping indel records merged): {n_events}")
    R(f"  log10 LR = {log10lr:.2f}   (>=3 strong, 1-3 moderate, <1 uninformative)")
    if n_priv_shared >= 2 or log10lr >= 3:
        verdict = "CLONALLY RELATED (same neoplasm / subclone / metastasis)"
    elif n_priv_shared == 1 or log10lr >= 1:
        verdict = "PROBABLY RELATED - confirm with WES / methylation / breakpoint matching"
    elif len(shared) == 0 and (len(p1) + len(p2)) >= 4:
        verdict = "NO SHARED SOMATIC VARIANTS - favours two independent primaries"
    else:
        verdict = "INDETERMINATE - too few informative variants; do WES or methylation array"
    R(f"  Verdict: {verdict}")
    if shared.rescued_in_T2.sum() > 0:
        R(f"  NOTE: {int(shared.rescued_in_T2.sum())} shared variant(s) were NOT called in T2 by the "
          f"pipeline but are present in T2 reads; the T2 caller was limited by purity {a.purity2}.")
    R()

    # ---- Step 7: 9p21 ------------------------------------------------------------
    R("-" * 78)
    R("STEP 7  9p21 (CDKN2A / CDKN2B / MTAP) COPY NUMBER")
    if a.targets:
        log("Step 7/8: 9p21 coverage")
        targets = read_bed(a.targets)
        res9 = analyse_9p21(bam1, bam2, targets, a.genome, a.purity1, a.purity2, a.out, a.min_mq)
        if res9 is not None:
            _, summ, bp = res9
            R("  gene     T1 ratio  T1 abs CN   T2 ratio  T2 abs CN   (purity-corrected)")
            for g in ("MTAP", "CDKN2A", "CDKN2B"):
                r1, c1 = summ[("T1", g)]
                r2, c2 = summ[("T2", g)]
                R(f"  {g:<8s} {r1:8.2f}  {c1:9.2f}   {r2:8.2f}  {c2:9.2f}")
            R(f"  deleted-target span T1: {bp['T1']}")
            R(f"  deleted-target span T2: {bp['T2']}")
            if bp["T1"] and bp["T2"]:
                same = abs(bp["T1"][0] - bp["T2"][0]) < 20000 and abs(bp["T1"][1] - bp["T2"][1]) < 20000
                R("  Deletion boundaries " + ("MATCH at target resolution -> supports common origin"
                                              if same else "DIFFER -> examine cn_9p21_profile.png; consider SNP array/WGS"))
            R("  Plot: cn_9p21_profile.png ; table: cn_9p21_targets.tsv")
    else:
        R("  skipped (provide --targets BED). Manual check with the ratios you already have:")
    # always give the purity reconciliation using user-supplied ratios if desired
    R("  Purity reconciliation formula: absCN = (2*ratio - 2*(1-purity)) / purity")
    R(f"    e.g. ratio 0.32 at purity {a.purity1:.2f} -> CN {cn_from_ratio(0.32, a.purity1):.2f}; "
      f"ratio 0.77 at purity {a.purity2:.2f} -> CN {cn_from_ratio(0.77, a.purity2):.2f}")
    R("    (both ~0 => same homozygous deletion seen through different purity)")
    R()

    # ---- Step 8: allelic imbalance --------------------------------------------
    R("-" * 78)
    R("STEP 8  ALLELIC IMBALANCE / LOH LANDSCAPE (germline het SNPs, mirrored BAF by arm)")
    if ident is not None and ident["arm_table"] is not None and not ident["arm_table"].empty:
        arm = ident["arm_table"]
        R(f"  arms evaluated: {len(arm)}   correlation of mBAF(T1) vs mBAF(T2): {ident['arm_corr']:.2f}")
        top = arm.sort_values("mbaf_T1", ascending=False).head(8)
        R("  arm    n   mBAF_T1  mBAF_T2")
        for _, r in top.iterrows():
            R(f"  {r.arm:<5s} {int(r.n):4d}   {r.mbaf_T1:.3f}    {r.mbaf_T2:.3f}")
        R("  Shared arms with high mBAF in both samples = shared LOH events (clonal evidence).")
        R("  Low mBAF in T2 across the board is expected at low purity and is NOT evidence against.")
        R("  Table: allelic_imbalance_by_arm.tsv")
    else:
        R("  not enough het SNPs.")
    R()

    # ---- Interpretation guide --------------------------------------------------
    R("-" * 78)
    R("HOW TO READ THIS")
    R("  * >=2 shared PRIVATE (non-hotspot) somatic variants with matching CCF = one neoplasm.")
    R("  * Shared hotspot(s) only = weak evidence; go to WES / methylation array / 9p21 breakpoint.")
    R("  * Zero shared variants with >=4 confident private variants in each = independent primaries,")
    R("    PROVIDED T2 purity is high enough that a clonal variant would have been seen")
    R(f"    (at purity {a.purity2:.2f} a clonal heterozygous variant is expected at VAF ~{a.purity2/2:.2f}).")
    R("  * GIS, TMB, MSI and CDKN2A ratio differences track purity, not biology; do not use them.")
    R("  * HLA type is germline; identical HLA is uninformative. HLA LOH (LOHHLA) can be subclonal.")
    R("  * Lineage call: SF-1 (ACC) vs WT1/D2-40/CK5-6/BAP1-loss (mesothelioma); shared CTNNB1/ZNRF3")
    R("    favours ACC, shared BAP1/NF2 favours mesothelioma. Methylation classifier settles lineage.")
    R("=" * 78)

    with open(os.path.join(a.out, "clonality_report.txt"), "w") as fh:
        fh.write("\n".join(report) + "\n")
    log(f"done. Outputs in {a.out}/")


if __name__ == "__main__":
    main()
