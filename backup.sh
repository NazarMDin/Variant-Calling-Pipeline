#!/bin/bash
# Manual Docker fallback script for variant calling and benchmarking
# Run from the project root directory (where FASTQ or BAM resides)
# Usage: bash backup_script.sh

set -e   # exit on error

# ---------- Configuration ----------
BASE="/home/kashif/nazar/acb"
REF="$BASE/GCF_000001405.40_GRCh38.p14_genomic.fna"
TRUTH_VCF="$BASE/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="$BASE/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BED="$BASE/chr1-22.bed"
OUTDIR="$BASE/results/manual"
mkdir -p "$OUTDIR"/{fastq,bam,deepvariant,clair3,renamed_vcfs,chr_reference,benchmark}

# Number of threads
THREADS=8

# Containers
SAMTOOLS="staphb/samtools:latest"
MINIMAP2="staphb/minimap2:latest"
PICARD="broadinstitute/picard:latest"
DEEPVARIANT="google/deepvariant:1.6.0"
CLAIR3="hkubal/clair3:latest"
BCFTOOLS="biocontainers/bcftools:v1.9-1-deb_cv1"
HAPPY="quay.io/biocontainers/hap.py:0.3.7--py27_0"

# ---------- Helper function ----------
run_docker() {
    local img=$1; shift
    docker run --rm -u $(id -u):$(id -g) -v "$BASE":"$BASE" -w "$BASE" "$img" "$@"
}

# ---------- Step 1: Ensure reference indexes ----------
if [ ! -f "$REF.fai" ]; then
    echo "Indexing reference..."
    run_docker $SAMTOOLS samtools faidx "$REF"
fi
dictFile="${REF%.fna}.dict"
if [ ! -f "$dictFile" ]; then
    echo "Creating reference dictionary with Picard..."
    run_docker $PICARD picard CreateSequenceDictionary R="$REF" O="$dictFile"
fi

# ---------- Step 2: Create chr1-22 BED ----------
echo "Generating autosome BED..."
grep -E '^NC_0000(0[1-9]|1[0-9]|2[0-2])\.' "$REF.fai" \
    | cut -f1,2 | awk -v OFS='\t' '{print $1, 0, $2}' > "$BED"

# ---------- Step 3: Alignment (if needed) ----------
if [ ! -f "$BASE/aligned_reads.bam" ]; then
    echo "Aligning reads with minimap2..."
    FASTQ="$BASE/m21009_241011_231051.hifi_reads.fastq.gz"   # adjust if different
    run_docker $MINIMAP2 sh -c "minimap2 -ax map-hifi -t $THREADS $REF $FASTQ | samtools sort -@ $THREADS -o $BASE/aligned_reads.bam -"
    run_docker $SAMTOOLS samtools index "$BASE/aligned_reads.bam"
else
    echo "Using existing aligned BAM: $BASE/aligned_reads.bam"
fi

BAM="$BASE/aligned_reads.bam"
BAI="$BAM.bai"

# ---------- Step 4: DeepVariant (using run_deepvariant) ----------
echo "Running DeepVariant..."
run_docker $DEEPVARIANT /opt/deepvariant/bin/run_deepvariant \
    --model_type PACBIO \
    --ref "$REF" \
    --reads "$BAM" \
    --regions "$BED" \
    --output_vcf "$OUTDIR/deepvariant/deepvariant.vcf.gz" \
    --num_shards $THREADS

# ---------- Step 5: Clair3 ----------
echo "Running Clair3..."
run_docker $CLAIR3 run_clair3.sh \
    --bam_fn="$BAM" \
    --ref_fn="$REF" \
    --threads=$THREADS \
    --platform="hifi" \
    --model_path="/opt/models/hifi" \
    --bed_fn="$BED" \
    --output="$OUTDIR/clair3"
cp "$OUTDIR/clair3/merge_output.vcf.gz" "$OUTDIR/clair3/clair3.vcf.gz"

# ---------- Step 6: Create contig mapping ----------
echo "Generating contig mapping..."
awk '$1 ~ /^NC_0000/ {split($1,a,"."); n=a[1]; sub(/NC_0000/,"",n); n=int(n); if(n>=1 && n<=22) print $1, "chr"n}' "$REF.fai" > "$BASE/accession_to_chr.map"

# ---------- Step 7: Rename VCFs ----------
echo "Renaming DeepVariant VCF..."
run_docker $BCFTOOLS bcftools annotate --rename-chrs "$BASE/accession_to_chr.map" -Oz -o "$OUTDIR/renamed_vcfs/deepvariant.renamed.vcf.gz" "$OUTDIR/deepvariant/deepvariant.vcf.gz"

echo "Renaming Clair3 VCF..."
run_docker $BCFTOOLS bcftools annotate --rename-chrs "$BASE/accession_to_chr.map" -Oz -o "$OUTDIR/renamed_vcfs/clair3.renamed.vcf.gz" "$OUTDIR/clair3/clair3.vcf.gz"

# ---------- Step 8: Create chr-prefixed reference ----------
echo "Creating chr-prefixed reference..."
awk 'NR==FNR {map[$1]=$2; next}
     /^>/ {
         split($0,a," ");
         h=substr(a[1],2);
         if(h in map) {
             print ">"map[h] " " substr($0, index($0,$2));
         } else {
             print $0
         }
         next
     }
     { print }' "$BASE/accession_to_chr.map" "$REF" > "$OUTDIR/chr_reference/chr_ref.fna"
run_docker $SAMTOOLS samtools faidx "$OUTDIR/chr_reference/chr_ref.fna"
run_docker $PICARD picard CreateSequenceDictionary R="$OUTDIR/chr_reference/chr_ref.fna" O="$OUTDIR/chr_reference/chr_ref.dict"

# ---------- Step 9: Benchmark ----------
echo "Benchmarking DeepVariant..."
run_docker $HAPPY /usr/local/bin/hap.py \
    "$TRUTH_VCF" \
    "$OUTDIR/renamed_vcfs/deepvariant.renamed.vcf.gz" \
    -f "$TRUTH_BED" \
    -r "$OUTDIR/chr_reference/chr_ref.fna" \
    -o "$OUTDIR/benchmark/happy_deepvariant" \
    --engine=vcfeval \
    --threads=$THREADS

echo "Benchmarking Clair3..."
run_docker $HAPPY /usr/local/bin/hap.py \
    "$TRUTH_VCF" \
    "$OUTDIR/renamed_vcfs/clair3.renamed.vcf.gz" \
    -f "$TRUTH_BED" \
    -r "$OUTDIR/chr_reference/chr_ref.fna" \
    -o "$OUTDIR/benchmark/happy_clair3" \
    --engine=vcfeval \
    --threads=$THREADS

echo "All done! Results in $OUTDIR/benchmark/"