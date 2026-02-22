minimap2 -ax map-hifi -t 8 \
  GCF_000001405.40_GRCh38.p14_genomic.fna \
  m21009_241011_231051.hifi_reads.fastq.gz \
  | samtools sort -@ 4 -o aligned_reads.bam -
samtools index aligned_reads.bam