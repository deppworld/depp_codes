#!/usr/bin/env bash
# run_clonality.sh -- end-to-end wrapper. Edit the VARIABLES block, then: bash run_clonality.sh
# Everything runs locally. No PHI leaves the machine.
set -euo pipefail

# ------------------------------- VARIABLES ---------------------------------- #
T1_BAM=T1.bam                       # TSO500 DRAGEN BAM (+ .bai next to it)
T2_BAM=T2.bam
T1_VCF=T1.hard-filtered.vcf.gz      # or *_MergedSmallVariants.genome.vcf.gz (gVCF also OK)
T2_VCF=T2.hard-filtered.vcf.gz
NORMAL_VCF=""                       # blood / normal VCF if you have one (strongly recommended)
T1_CVO=T1_CombinedVariantOutput.tsv # optional; gives gene / p. annotation
T2_CVO=T2_CombinedVariantOutput.tsv
TARGETS=TSO500_targets.bed          # panel manifest BED (hg19 for standard TSO500)
GENOME=hg19
PURITY1=0.75                        # pathologist estimate or ploidy-tool estimate
PURITY2=0.30
OUT=clonality_out
THREADS=8
# ---------------------------------------------------------------------------- #

# 0. Optional: if you only have FASTQs, align them (needs bwa + samtools + ref fasta).
#    TSO500 BAMs from DRAGEN are already fine; skip this block.
# REF=hs37d5.fa
# for S in T1 T2; do
#   bwa mem -t $THREADS -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA" $REF ${S}_R1.fastq.gz ${S}_R2.fastq.gz \
#     | samtools sort -@ $THREADS -o $S.bam - && samtools index $S.bam
# done

# 1. Sanity: indexes present, contig naming consistent
for B in "$T1_BAM" "$T2_BAM"; do
  [ -f "$B.bai" ] || [ -f "${B%.bam}.bai" ] || samtools index "$B"
done
for V in "$T1_VCF" "$T2_VCF"; do
  [ -f "$V.tbi" ] || { [[ "$V" == *.gz ]] && tabix -p vcf "$V"; } || true
done
echo "BAM contigs T1: $(samtools view -H "$T1_BAM" | grep -c '^@SQ') ; T2: $(samtools view -H "$T2_BAM" | grep -c '^@SQ')"

# 2. Main comparison
ARGS=(--t1-vcf "$T1_VCF" --t2-vcf "$T2_VCF" --t1-bam "$T1_BAM" --t2-bam "$T2_BAM"
      --purity1 "$PURITY1" --purity2 "$PURITY2" --genome "$GENOME" --out "$OUT")
[ -n "$NORMAL_VCF" ] && ARGS+=(--normal-vcf "$NORMAL_VCF")
[ -f "$T1_CVO" ]     && ARGS+=(--t1-cvo "$T1_CVO")
[ -f "$T2_CVO" ]     && ARGS+=(--t2-cvo "$T2_CVO")
[ -f "$TARGETS" ]    && ARGS+=(--targets "$TARGETS")
python3 "$(dirname "$0")/clonality_check.py" "${ARGS[@]}"

# 3. Optional: formal statistic with the Clonality R package (Ostrovnaya et al.)
#    Rscript -e 'install.packages("BiocManager"); BiocManager::install("Clonality")'
cat > "$OUT/clonality_R_input.R" <<'EOF'
# Builds a 0/1 mutation matrix from variants_T1_T2.tsv and runs the LR test.
library(Clonality)
v <- read.delim("variants_T1_T2.tsv")
v <- v[v$status %in% c("SHARED","T1_PRIVATE","T2_PRIVATE") & !v$chip_flag, ]
m <- rbind(T1 = as.integer(v$status %in% c("SHARED","T1_PRIVATE")),
           T2 = as.integer(v$status %in% c("SHARED","T2_PRIVATE")))
colnames(m) <- v$key
# freq = per-variant probability of being mutated in an unrelated tumor.
# Use hotspot recurrence for hotspots, and a tiny value for private variants.
freq <- ifelse(is.na(v$hotspot_freq), 1e-6, v$hotspot_freq)
res <- mutation.proba(m, freq)   # returns P(same clone) style probabilities
print(res)
EOF
echo "Formal test (optional): cd $OUT && Rscript clonality_R_input.R"
echo "Report: $OUT/clonality_report.txt"
