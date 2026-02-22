#!/bin/bash
# Manual Docker commands to run variant calling and benchmarking
# Run this from the project root directory

# Set paths (adjust if needed)
BAM="aligned_reads.bam"
REF="GCF_000001405.40_GRCh38.p14_genomic.fna"
BED="chr1-22.bed"
OUTDIR_BASE="results/manual"
mkdir -p $OUTDIR_BASE/deepvariant $OUTDIR_BASE/clair3 $OUTDIR_BASE/renamed_vcfs $OUTDIR_BASE/benchmark

# DeepVariant: make_examples
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/make_examples \
  --mode calling \
  --ref /data/$REF \
  --reads /data/$BAM \
  --examples /data/$OUTDIR_BASE/deepvariant/examples.tfrecord.gz \
  --regions /data/$BED \
  --task 0

# DeepVariant: call_variants
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/call_variants \
  --outfile /data/$OUTDIR_BASE/deepvariant/call_variants_output.tfrecord.gz \
  --examples /data/$OUTDIR_BASE/deepvariant/examples.tfrecord.gz \
  --checkpoint /opt/models/wgs

# DeepVariant: postprocess_variants
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/postprocess_variants \
  --ref /data/$REF \
  --infile /data/$OUTDIR_BASE/deepvariant/call_variants_output.tfrecord.gz \
  --outfile /data/$OUTDIR_BASE/deepvariant/deepvariant.vcf.gz

# Clair3
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  hkubal/clair3:latest \
  run_clair3.sh \
  --bam_fn=/data/$BAM \
  --ref_fn=/data/$REF \
  --threads=8 \
  --platform="hifi" \
  --model_path=/opt/models/hifi \
  --bed_fn=/data/$BED \
  --output=/data/$OUTDIR_BASE/clair3

cp $OUTDIR_BASE/clair3/merge_output.vcf.gz $OUTDIR_BASE/clair3/clair3.vcf.gz

# Generate mapping file
awk '$1 ~ /^NC_0000/ {split($1,a,"."); n=a[1]; sub(/NC_0000/,"",n); n=int(n); if(n>=1 && n<=22) print $1, "chr"n}' $REF.fai > accession_to_chr.map

# Rename contigs
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools annotate --rename-chrs /data/accession_to_chr.map -Oz -o /data/$OUTDIR_BASE/renamed_vcfs/deepvariant.renamed.vcf.gz /data/$OUTDIR_BASE/deepvariant/deepvariant.vcf.gz

docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  biocontainers/bcftools:v1.9-1-deb_cv1 \
  bcftools annotate --rename-chrs /data/accession_to_chr.map -Oz -o /data/$OUTDIR_BASE/renamed_vcfs/clair3.renamed.vcf.gz /data/$OUTDIR_BASE/clair3/clair3.vcf.gz

# Benchmark with hap.py (using quay.io image)
HAPPY_IMAGE="quay.io/biocontainers/hap.py:0.3.7--py27_0"
TRUTH_VCF="/data/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="/data/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# DeepVariant
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  $HAPPY_IMAGE \
  /usr/local/bin/hap.py \
  $TRUTH_VCF \
  /data/$OUTDIR_BASE/renamed_vcfs/deepvariant.renamed.vcf.gz \
  -f $TRUTH_BED \
  -r /data/$REF \
  -o /data/$OUTDIR_BASE/benchmark/happy_deepvariant \
  --engine=vcfeval \
  --threads=8

# Clair3
docker run --rm -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  $HAPPY_IMAGE \
  /usr/local/bin/hap.py \
  $TRUTH_VCF \
  /data/$OUTDIR_BASE/renamed_vcfs/clair3.renamed.vcf.gz \
  -f $TRUTH_BED \
  -r /data/$REF \
  -o /data/$OUTDIR_BASE/benchmark/happy_clair3 \
  --engine=vcfeval \
  --threads=8

echo "All done! Results in $OUTDIR_BASE/benchmark/"