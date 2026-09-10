#!/usr/bin/env python3
"""Build a small synthetic T1/T2 dataset to smoke-test clonality_check.py.
Scenario: same clone, T1 purity 0.75, T2 purity 0.30; CDKN2A homozygous deletion in both;
T2 caller misses a low-VAF shared variant (to be rescued by force-calling)."""
import os, random, sys
import pysam
random.seed(7)
out = sys.argv[1] if len(sys.argv) > 1 else "synthetic"
os.makedirs(out, exist_ok=True)

CONTIGS = [("1", 250_000_000), ("3", 200_000_000), ("9", 141_000_000), ("17", 81_000_000)]
BASES = "ACGT"

def randseq(n):
    return "".join(random.choice(BASES) for _ in range(n))

# reference context per site (we just need consistent ref bases)
sites = {  # key -> (chrom,pos,ref,alt, vaf_T1, vaf_T2, gene, pdot, kind)
    "s1": ("17", 7_577_120, "C", "T", 0.60, 0.20, "TP53",   "p.(R273H)", "hotspot"),
    "s2": ("3",  41_266_137, "C", "T", 0.35, 0.14, "CTNNB1", "p.(S45F)",  "hotspot"),
    "s3": ("1",  115_256_529, "T", "C", 0.36, 0.15, "NRAS",  "p.(Q61R)",  "hotspot"),
    "s4": ("1",  27_100_000, "G", "A", 0.37, 0.13, "ARID1A", "p.(E1000K)", "private"),
    "s5": ("9",  35_000_000, "A", "G", 0.34, 0.04, "MLLT3",  "p.(P45L)",  "private_missed_in_T2"),
    "s6": ("17", 41_244_000, "T", "G", 0.20, 0.00, "BRCA1",  "p.(K1000N)", "T1_private_subclonal"),
    "s7": ("3",  178_936_091, "G", "A", 0.00, 0.10, "PIK3CA", "p.(E545K)", "T2_private"),
    "g1": ("1",  100_000_000, "A", "C", 0.50, 0.50, "RPL5",  "p.(=)",     "germline"),
}
# germline het SNPs for identity/BAF (same genotype in both)
het = [(c, 5_000_000 + i * 2_000_000, "A", "G") for c, _ in CONTIGS for i in range(6)]

def write_bam(path, vaf_idx, depth=300, purity=0.75, cdkn2a_ratio=0.32):
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": c, "LN": ln} for c, ln in CONTIGS]}
    reads = []
    def add_reads(chrom, pos, ref, alt, vaf, n):
        for i in range(n):
            start = pos - 1 - random.randint(20, 80)
            seq = list(randseq(101))
            off = pos - 1 - start
            seq[off] = alt if random.random() < vaf else ref
            a = pysam.AlignedSegment()
            a.query_name = f"r_{chrom}_{pos}_{i}"
            a.query_sequence = "".join(seq)
            a.flag = 0 if i % 2 else 16
            a.reference_id = [c for c, _ in CONTIGS].index(chrom)
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = ((0, 101),)
            a.query_qualities = pysam.qualitystring_to_array("I" * 101)
            reads.append(a)
    for k, (c, p, r, alt, v1, v2, *_ ) in sites.items():
        add_reads(c, p, r, alt, (v1, v2)[vaf_idx], depth)
    for c, p, r, alt in het:
        add_reads(c, p, r, alt, 0.5, 120)
    # 9p21 targets: normal targets and deleted targets
    for c, p, r, alt in targets_norm:
        add_reads(c, p, r, alt, 0.0, depth)
    for c, p, r, alt in targets_del:
        add_reads(c, p, r, alt, 0.0, int(depth * cdkn2a_ratio))
    reads.sort(key=lambda a: (a.reference_id, a.reference_start))
    with pysam.AlignmentFile(path, "wb", header=hdr) as fh:
        for a in reads:
            fh.write(a)
    pysam.index(path)

# target BED: 30 off-chr9 targets + 9p21 targets (MTAP, CDKN2A, CDKN2B, flanks)
targets_norm = [(c, 10_500_000 + i * 3_000_000, "A", "G") for c in ("1", "3", "17") for i in range(10)]
p21 = [21_750_000, 21_830_000, 21_850_000, 21_975_000, 21_985_000, 21_993_000,
       22_004_000, 22_007_000, 22_100_000, 22_250_000]
deleted = {21_830_000, 21_850_000, 21_975_000, 21_985_000, 21_993_000, 22_004_000, 22_007_000}
targets_del = [("9", p, "A", "G") for p in p21 if p in deleted]
targets_norm += [("9", p, "A", "G") for p in p21 if p not in deleted]
with open(f"{out}/targets.bed", "w") as fh:
    for c, p, *_ in targets_norm + targets_del:
        fh.write(f"chr{c}\t{p-60}\t{p+60}\ttarget_{c}_{p}\n")

write_bam(f"{out}/T1.bam", 0, purity=0.75, cdkn2a_ratio=0.32)
write_bam(f"{out}/T2.bam", 1, purity=0.30, cdkn2a_ratio=0.77)

def write_vcf(path, idx, drop_keys):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, ln in CONTIGS:
            fh.write(f"##contig=<ID=chr{c},length={ln}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n')
        fh.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="">\n')
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="">\n')
        fh.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        rows = []
        for k, (c, p, r, alt, v1, v2, *_ ) in sites.items():
            v = (v1, v2)[idx]
            if v == 0 or k in drop_keys:
                continue
            dp = 300; ad = int(dp * v)
            rows.append((c, p, r, alt, dp, ad, v))
        for c, p, r, alt in het:
            rows.append((c, p, r, alt, 120, 60, 0.5))
        rows.sort(key=lambda x: ([c for c, _ in CONTIGS].index(x[0]), x[1]))
        for c, p, r, alt, dp, ad, v in rows:
            fh.write(f"chr{c}\t{p}\t.\t{r}\t{alt}\t100\tPASS\t.\tGT:AD:DP:AF\t0/1:{dp-ad},{ad}:{dp}:{v:.4f}\n")
    pysam.tabix_compress(path, path + ".gz", force=True)
    pysam.tabix_index(path + ".gz", preset="vcf", force=True)

write_vcf(f"{out}/T1.vcf", 0, set())
write_vcf(f"{out}/T2.vcf", 1, {"s5"})   # caller missed s5 in T2

def write_cvo(path):
    with open(path, "w") as fh:
        fh.write("[Analysis Details]\nfoo\tbar\n\n[Small Variants]\n")
        fh.write("Gene\tChromosome\tGenomic Position\tReference Call\tAlternative Call\tAllele Frequency\tDepth\tP-Dot Notation\tC-Dot Notation\tConsequence(s)\tAffected Exon(s)\n")
        for k, (c, p, r, alt, v1, v2, gene, pdot, kind) in sites.items():
            fh.write(f"{gene}\tchr{c}\t{p}\t{r}\t{alt}\t{v1}\t300\t{pdot}\tc.?\tmissense_variant\t.\n")
        fh.write("\n[Copy Number Variants]\n")
write_cvo(f"{out}/T1_CombinedVariantOutput.tsv")
write_cvo(f"{out}/T2_CombinedVariantOutput.tsv")
print("synthetic data written to", out)
