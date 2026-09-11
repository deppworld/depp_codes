#!/usr/bin/env python3
"""
estimate_purity.py -- estimate tumor purity for a TSO500 sample from its own data,
for use when no pathology tumor-cell percentage is available.

Three independent estimators are reported; use the pathology number if you later
get one, otherwise the consensus printed at the bottom.

  A. Homozygous-deletion method : for a gene at copy number 0 in tumor cells, the
     coverage fold-change equals (1 - purity).  purity = 1 - ratio.
     Uses --del-ratio values you pass (e.g. CDKN2A 0.32) or reads them from the
     DRAGEN CNV VCF (--cnv-vcf) for --del-genes.
  B. Somatic VAF method         : a clonal heterozygous mutation in a diploid region
     sits at VAF = purity/2.  Somatic-looking variants (not in 0.40-0.60 or >0.90
     germline bands, PASS, depth >= 50) are clustered; purity = 2 x the upper
     cluster of VAFs (the clonal peak).  Reported as median of top cluster and as
     2 x 90th percentile (upper bound).
  C. TP53 / LOH method          : a clonal mutation with loss of the wild-type allele
     sits at VAF = purity.  If a TP53 (or other named gene) mutation is present,
     its VAF is a direct purity estimate (upper bound if no LOH).

Usage:
  python3 estimate_purity.py --vcf T1.dna.hard-filtered.vcf \
        --cvo T1_CombinedVariantOutput.tsv \
        --del-ratio CDKN2A=0.32 --del-ratio CDKN2B=0.36
  # or let it pull fold-changes from the CNV VCF:
  python3 estimate_purity.py --vcf T1.dna.hard-filtered.vcf --cnv-vcf T1.dna.cnv.vcf
"""
import argparse
import re
import sys

import numpy as np
import pysam

sys.path.insert(0, __import__("os").path.dirname(__file__))
try:
    from clonality_check import parse_cvo, norm_chrom  # reuse the CVO parser
except Exception:  # standalone fallback
    parse_cvo = None

    def norm_chrom(c):
        c = str(c)
        return c[3:] if c.lower().startswith("chr") else c


def somatic_vafs(vcf_path, min_dp=50, min_alt=5):
    vf = pysam.VariantFile(vcf_path)
    out = []
    for rec in vf:
        if rec.alts is None or rec.alts[0] in (None, "<NON_REF>"):
            continue
        filt = list(rec.filter.keys())
        if filt and filt != ["PASS"]:
            continue
        if not len(vf.header.samples):
            continue
        s = rec.samples[0]
        dp = s.get("DP")
        af = s.get("AF")
        ad = s.get("AD")
        if af is None and ad and dp:
            af = ad[1] / dp
        if af is None or dp is None:
            continue
        af = af[0] if isinstance(af, tuple) else af
        if dp < min_dp:
            continue
        if ad and len(ad) > 1 and ad[1] < min_alt:
            continue
        # drop germline-looking VAFs
        if 0.40 <= af <= 0.60 or af > 0.90:
            continue
        out.append((f"{norm_chrom(rec.chrom)}:{rec.pos}:{rec.ref}:{rec.alts[0]}", float(af), int(dp)))
    return out


def cnv_fold_changes(cnv_vcf, genes):
    """Pull per-gene fold change from a DRAGEN TSO500 CNV VCF (INFO FC / FOLD_CHANGE or SM)."""
    res = {}
    vf = pysam.VariantFile(cnv_vcf)
    for rec in vf:
        info = dict(rec.info)
        gene = None
        for k in ("GENE", "Gene", "ANN", "SVTYPE"):
            pass
        # DRAGEN TSO500 writes the gene in ID or INFO; be permissive
        text = f"{rec.id} {info}"
        for g in genes:
            if re.search(rf"\b{g}\b", text):
                gene = g
        if gene is None:
            continue
        fc = None
        for k in ("FC", "FOLD_CHANGE", "FoldChange", "SM"):
            if k in info:
                fc = info[k]
                fc = fc[0] if isinstance(fc, tuple) else fc
                break
        if fc is not None:
            res[gene] = float(fc)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="hard-filtered small-variant VCF")
    ap.add_argument("--cvo", default=None, help="CombinedVariantOutput.tsv (for gene names)")
    ap.add_argument("--cnv-vcf", default=None, help="DRAGEN CNV VCF to read fold changes from")
    ap.add_argument("--del-ratio", action="append", default=[],
                    help="GENE=ratio for a homozygously deleted gene, e.g. CDKN2A=0.32 (repeatable)")
    ap.add_argument("--del-genes", default="CDKN2A,CDKN2B,MTAP,PTEN,RB1,SMAD4",
                    help="genes to look up in --cnv-vcf")
    ap.add_argument("--loh-genes", default="TP53", help="genes whose mutation VAF ~ purity when LOH present")
    a = ap.parse_args()

    print("=" * 70)
    print(f"PURITY ESTIMATE  {a.vcf}")
    print("=" * 70)

    # ---- A. deletion method --------------------------------------------------
    ratios = {}
    for d in a.del_ratio:
        g, r = d.split("=")
        ratios[g.strip()] = float(r)
    if a.cnv_vcf:
        ratios.update(cnv_fold_changes(a.cnv_vcf, a.del_genes.split(",")))
    estA = []
    print("\nA. Homozygous-deletion method (purity = 1 - fold change)")
    if not ratios:
        print("   no deletion ratios given (use --del-ratio CDKN2A=0.32 or --cnv-vcf)")
    for g, r in ratios.items():
        if r < 0.85:
            p = 1 - r
            estA.append(p)
            print(f"   {g:<8s} ratio {r:.2f}  ->  purity {p:.2f}   (valid only if deletion is homozygous & clonal)")
        else:
            print(f"   {g:<8s} ratio {r:.2f}  ->  not deleted enough to use")

    # ---- B. VAF method ---------------------------------------------------------
    vafs = somatic_vafs(a.vcf)
    ann = {}
    if a.cvo and parse_cvo:
        c = parse_cvo(a.cvo)
        ann = {r.key: (r.gene, r.pdot) for _, r in c.iterrows()}
    print(f"\nB. Somatic VAF method ({len(vafs)} somatic-looking PASS variants, depth>=50)")
    estB = None
    if len(vafs) >= 1:
        v = np.array(sorted([x[1] for x in vafs], reverse=True))
        # clonal peak = variants within 0.08 of the maximum (excluding an isolated outlier
        # if it is > 1.6x the next value, which suggests LOH/amplification)
        top = v[0]
        if len(v) > 1 and v[0] > 1.6 * v[1]:
            print(f"   highest VAF {v[0]:.3f} is an outlier vs next {v[1]:.3f} -> likely LOH/CN gain; "
                  f"treated under method C")
            top = v[1]
        peak = v[(v >= top - 0.08) & (v <= top + 0.001)]
        estB = 2 * float(np.median(peak))
        print(f"   clonal VAF peak (n={len(peak)}): median {np.median(peak):.3f}  ->  purity ~ {min(estB,1):.2f}")
        print(f"   2 x 90th percentile VAF        :  {min(2*np.percentile(v,90),1):.2f}  (upper bound)")
        print("   top variants:")
        for k, af, dp in sorted(vafs, key=lambda x: -x[1])[:8]:
            g, p = ann.get(k, ("", ""))
            print(f"     {k:<30s} VAF {af:.3f} DP {dp:5d}  {g} {p}")
    else:
        print("   no usable somatic variants")

    # ---- C. LOH gene method -----------------------------------------------------
    print(f"\nC. LOH-gene method (VAF of {a.loh_genes} mutation ~ purity if wild-type allele lost)")
    estC = []
    for k, af, dp in vafs:
        g, p = ann.get(k, ("", ""))
        if g in a.loh_genes.split(","):
            estC.append(af)
            print(f"   {g} {p}  VAF {af:.3f}  ->  purity ~ {af:.2f} (if LOH) or ~ {min(2*af,1):.2f} (if no LOH)")
    if not estC:
        print("   no mutation in these genes found (needs --cvo for gene names)")

    # ---- consensus -----------------------------------------------------------------
    print("\n" + "-" * 70)
    cands = [p for p in estA] + ([estB] if estB else [])
    if cands:
        cons = float(np.median(cands))
        lo, hi = min(cands), max(cands)
        print(f"CONSENSUS purity ~ {cons:.2f}   (range {lo:.2f}-{hi:.2f})")
        print(f"Suggested for clonality_check.py: --purityN {cons:.2f}   and re-run with {lo:.2f} and {hi:.2f} as sensitivity check")
        if hi - lo > 0.2:
            print("NOTE: estimators disagree by >0.20. Prefer the deletion-based value if the 9p21 deletion is")
            print("      clearly homozygous; prefer the VAF value if the top VAF cluster has >=3 variants.")
    else:
        print("No estimate possible; obtain pathology tumor-cell % from the macrodissected H&E.")
    print("-" * 70)


if __name__ == "__main__":
    main()
