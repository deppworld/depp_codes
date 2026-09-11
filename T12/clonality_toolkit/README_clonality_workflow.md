# Clonality workflow: T1 vs T2 (same neoplasm or two primaries?)

Everything below runs on your own workstation. No PHI is needed by the scripts beyond the
files themselves, and nothing is transmitted.

## Files

| file | purpose |
|---|---|
| `clonality_check.py` | the analysis (Python 3, needs `pysam pandas numpy matplotlib`) |
| `run_clonality.sh` | wrapper: edit the variables block, run once |
| `make_synthetic_test.py` | builds a toy T1/T2 dataset so you can verify the install before touching patient data |
 'estimate_purity.py'  estimate purity from the sequencing data
 	
Install: `pip install pysam pandas numpy matplotlib scipy` (samtools/tabix recommended).


Smoke test: `python3 make_synthetic_test.py demo && python3 clonality_check.py --t1-vcf demo/T1.vcf.gz --t2-vcf demo/T2.vcf.gz --t1-bam demo/T1.bam --t2-bam demo/T2.bam --purity1 0.75 --purity2 0.30 --targets demo/targets.bed --t1-cvo demo/T1_CombinedVariantOutput.tsv --t2-cvo demo/T2_CombinedVariantOutput.tsv --out demo/out`


## Estimate purity
python3 estimate_purity.py --vcf T1.dna.hard-filtered.vcf --cvo T1_CombinedVariantOutput.tsv --tmb-trace T1_TMB_Trace.tsv --del-ratio CDKN2A=0.32 --del-ratio CDKN2B=0.36

python3 estimate_purity.py --vcf T2.dna.hard-filtered.vcf --cvo T2_CombinedVariantOutput.tsv --tmb-trace T2_TMB_Trace.tsv --del-ratio CDKN2A=0.77 --del-ratio CDKN2B=0.84
## Run Script: 
Based on estimated purity (T1 70 and T2 25)

python3 clonality_check.py --t1-vcf T1.dna.hard-filtered.vcf --t2-vcf T2.dna.hard-filtered.vcf --t1-bam T1.dna.bam --t2-bam T2.dna.bam --t1-cvo T1_CombinedVariantOutput.tsv --t2-cvo T2_CombinedVariantOutput.tsv --t1-tmb-trace T1_TMB_Trace.tsv --t2-tmb-trace T2_TMB_Trace.tsv --targets TSO500_targets.bed --purity1 0.66 --purity2 0.20 --out results



## Inputs from the TSO500 run folder

* BAM + BAI for T1 and T2 (`Logs_Intermediates/StitchedRealigned/<sample>/<sample>.bam` or the DRAGEN `*.bam`)
* small-variant VCF per sample (`*.hard-filtered.vcf.gz` or `*_MergedSmallVariants.genome.vcf.gz`; the gVCF is fine and gives more heterozygous SNPs for the identity check)
* `<sample>_CombinedVariantOutput.tsv` (optional; supplies gene and p. annotation)
* the panel target BED (TSO500 manifest, hg19 unless you run a hg38 build)
* tumor purity for each sample: pathologist estimate from the H&E used for macrodissection, or from a ploidy tool. This is the single most important parameter, because most of the differences you are seeing (GIS 26 vs 2, CDKN2A ratio 0.32 vs 0.77) are purity effects.
* a germline VCF (blood or uninvolved tissue) if at all possible. Without it germline is inferred from VAF shape, which is weaker.

## Step-by-step

**Step 0 — Same patient?** The script force-genotypes heterozygous SNPs in both BAMs. Hom-ref vs hom-alt discordances should be ~0. If they are not, stop: sample swap or two patients. Do this before anything else; it is the cheapest way to avoid a catastrophic report.

**Step 1 — Harmonise the variant lists.** Both VCFs are parsed to `chrom:pos:ref:alt` keys. PASS-only by default (`--no-pass-only` to relax). Depth ≥ 20, alt reads ≥ 3, VAF ≥ 2 % for a caller call to count.

**Step 2 — Remove what is not somatic.** Variants in the normal VCF, variants sitting at 40–60 % (or > 90 %) in *both* samples without a normal, and low-VAF calls in CHIP genes with identical VAF in both samples are flagged. Germline-like calls are excluded from the comparison. CHIP flags are excluded but listed so you can review them.

**Step 3 — Force-call everything in both BAMs.** This is the part that matters for your case. T2 has low purity, so the DRAGEN caller may have dropped variants that are actually present at 2–5 % VAF. Every variant seen in either sample is re-counted by pileup in both BAMs (base quality ≥ 20, mapping quality ≥ 20). A variant is "present" with ≥ 3 alt reads and VAF ≥ 0.5 % (`--force-min-alt`, `--force-min-vaf`). Variants rescued this way are marked `rescued_in_T2` in the table.

**Step 4 — Classify.** SHARED / T1_PRIVATE / T2_PRIVATE, and HOTSPOT vs PRIVATE (non-recurrent). The built-in hotspot list covers the common ones (TP53 R175/R248/R273 etc., KRAS/NRAS G12/G13/Q61, BRAF V600, PIK3CA E542/E545/H1047, CTNNB1 D32–S45, IDH1/2, TERT promoter). Supply `--hotspots` to extend it.

**Step 5 — Cancer-cell fraction.** `CCF = VAF × (purity × CN + 2(1 − purity)) / purity`, CN = 2 unless you pass `--cn gene_cn.tsv` (columns `gene cn_T1 cn_T2`, taken from the TSO500 CNV output). CCF ≥ 0.7 in both = truncal. Shared truncal variants plus divergent private variants is exactly what a subclone or metastasis looks like.

**Step 6 — Likelihood ratio.** For each shared variant, LR contribution = 1 / P(chance match). Private variants get P ≈ (somatic burden / panel size) / 3, hotspots get their recurrence frequency. log10 LR ≥ 3 is strong support for a common origin. This is a simplified version of the Ostrovnaya method; `run_clonality.sh` writes an R script for the Bioconductor `Clonality` package if you want a citable statistic in the report.

**Step 7 — 9p21 deletion breakpoints.** Per-target median depth across MTAP–CDKN2A–CDKN2B, normalised to off-chr9 targets, for both BAMs, plus purity-corrected absolute copy number: `absCN = (2·ratio − 2(1 − purity)) / purity`. With your numbers: ratio 0.32 at ~75 % purity → CN ≈ 0.2; ratio 0.77 at ~30 % purity → CN ≈ 0.5. Both are a homozygous deletion seen through different purity. If the deleted-target span is identical in both samples that is a second, independent line of clonal evidence. Output: `cn_9p21_profile.png`, `cn_9p21_targets.tsv`. (Target-level resolution only; for base-pair breakpoints you need WGS or a SNP array.)

**Step 8 — LOH landscape.** Mirrored BAF of germline heterozygous SNPs per chromosome arm in both samples. Shared high-mBAF arms are shared LOH events. Low mBAF everywhere in T2 is expected at 30 % purity and is not evidence against clonality. This is also where HLA LOH would show up (chr6p); run LOHHLA if you need it formally.

## Decision rules for the report

1. ≥ 2 shared **private** somatic variants (or ≥ 1 private + matching 9p21 breakpoints) → one clonally related neoplasm. Report a single integrated diagnosis; describe T2 as a morphologically/immunophenotypically divergent subclone or metastatic deposit; state explicitly that GIS, TMB, MSI and CDKN2A-ratio differences are attributable to tumor purity.
2. Shared hotspot(s) only, no shared private variants → suggestive but not sufficient. Add WES (or the larger 500-gene comparison at gVCF level), or a DNA-methylation array, before committing.
3. No shared variants, ≥ 4 confident private variants in each, **and** T2 purity high enough that a clonal variant would have appeared (expected VAF ≈ purity/2) → two independent primaries.
4. Anything else → indeterminate; escalate to WES / methylation.

## Settling lineage once clonality is settled

* SF-1 is the discriminating marker: positive → adrenocortical; essentially never in mesothelioma. Calretinin is positive in a large fraction of ACC, so it must not carry a mesothelioma call on its own. Mesothelioma side: WT1, D2-40, CK5/6, strong diffuse keratin, BAP1 loss. MTAP-loss IHC will be lost in *both* given 9p21 deletion and is not discriminating here.
* Gene content of the shared variants: CTNNB1 / ZNRF3 / TERT / PRKAR1A / MEN1 → ACC; BAP1 / NF2 / SETD2 / LATS2 → mesothelioma.
* A methylation-array tissue-of-origin classifier (or RNA-seq) separates ACC from mesothelioma in one experiment and will also show whether T1 and T2 co-cluster.
* Clinically: ACC seeds pleura and peritoneum; a serosal metastasis with reactive mesothelial hyperplasia around it is a well-known mimic and would explain both the mesothelial immunophenotype and the low tumor content of T2.

## Things that are *not* evidence either way

Identical HLA type (germline), TMB 4.7 vs 3.9, MSI 4.9 % vs 6.7 % (both MSS, within noise), GIS 26 vs 2 and CDKN2A ratio 0.32 vs 0.77 (purity).

## Outputs

`clonality_report.txt` (readable summary), `variants_T1_T2.tsv` (every variant, both samples, force-called counts, CCF, flags), `shared_variants_LR.tsv`, `cn_9p21_targets.tsv`, `cn_9p21_profile.png`, `het_snps_T1_T2.tsv`, `allelic_imbalance_by_arm.tsv`, `clonality_R_input.R`.
