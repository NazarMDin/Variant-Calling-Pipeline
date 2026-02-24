#!/usr/bin/env nextflow

params.bam = "/home/kashif/nazar/acb/aligned_reads.bam"
params.ref = "/home/kashif/nazar/acb/GCF_000001405.40_GRCh38.p14_genomic.fna"
params.truth_vcf = "/home/kashif/nazar/acb/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truth_bed = "/home/kashif/nazar/acb/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
params.chr_bed = "/home/kashif/nazar/acb/chr1-22.bed"
params.outdir = "results"

params.singularity_deepvariant = "/home/kashif/nazar/acb/singularity_images/deepvariant_1.6.0.sif"
params.singularity_clair3 = "/home/kashif/nazar/acb/singularity_images/clair3_latest.sif"
params.singularity_happy = "/home/kashif/nazar/acb/singularity_images/happy.sif"

// Containers for auxiliary tools
params.bcftools_container = "biocontainers/bcftools:v1.9-1-deb_cv1"
params.samtools_container = "staphb/samtools:latest"

// Threads
params.num_shards = 8
params.clair3_threads = 8
params.happy_threads = 8

process deepvariant {
    container params.singularity_deepvariant
    publishDir "${params.outdir}/deepvariant", mode: 'copy'

    input:
    path bam
    path bai
    path ref
    path ref_fai
    path ref_dict
    path bed

    output:
    path "deepvariant.vcf.gz"

    script:
    """
    /opt/deepvariant/bin/make_examples \
      --mode calling \
      --ref "${ref}" \
      --reads "${bam}" \
      --examples "examples.tfrecord.gz" \
      --regions "${bed}" \
      --task 0

    /opt/deepvariant/bin/call_variants \
      --outfile "call_variants_output.tfrecord.gz" \
      --examples "examples.tfrecord.gz" \
      --checkpoint "/opt/models/wgs"

    /opt/deepvariant/bin/postprocess_variants \
      --ref "${ref}" \
      --infile "call_variants_output.tfrecord.gz" \
      --outfile "deepvariant.vcf.gz"
    """
}

process clair3 {
    container params.singularity_clair3
    publishDir "${params.outdir}/clair3", mode: 'copy'

    input:
    path bam
    path bai
    path ref
    path ref_fai
    path bed

    output:
    path "clair3.vcf.gz"

    script:
    """
    mkdir -p clair3_tmp
    run_clair3.sh \
      --bam_fn="${bam}" \
      --ref_fn="${ref}" \
      --threads=${params.clair3_threads} \
      --platform="hifi" \
      --model_path="/opt/models/hifi" \
      --bed_fn="${bed}" \
      --output="clair3_tmp"

    cp clair3_tmp/merge_output.vcf.gz clair3.vcf.gz
    """
}

process renameContigs {
    container params.bcftools_container
    publishDir "${params.outdir}/renamed_vcfs", mode: 'copy', pattern: "*.vcf.gz"

    input:
    tuple val(caller), path(vcf)
    path mapping_file

    output:
    tuple val(caller), path("renamed_${caller}.vcf.gz")

    script:
    """
    bcftools annotate --rename-chrs ${mapping_file} -Oz -o renamed_${caller}.vcf.gz ${vcf}
    bcftools index renamed_${caller}.vcf.gz
    """
}

process prepareChrRef {
    container params.samtools_container
    publishDir "${params.outdir}/chr_reference", mode: 'copy'

    input:
    path ref
    path mapping_file

    output:
    path "chr_ref.fna"
    path "chr_ref.fna.fai"
    path "chr_ref.dict"

    script:
    """
    awk 'NR==FNR {map[\$1]=\$2; next} /^>/ {h=\$1; sub(/^>/,"",h); if(h in map) {print ">"map[h]} else {print \$0}} !/^>/ {print}' \\
        ${mapping_file} ${ref} > chr_ref.fna
    samtools faidx chr_ref.fna
    samtools dict chr_ref.fna -o chr_ref.dict
    """
}

process hapPy {
    container params.singularity_happy
    publishDir "${params.outdir}/benchmark/${caller}", mode: 'copy'

    input:
    tuple val(caller), path(query_vcf)
    path truth_vcf
    path truth_bed
    path chr_ref
    path chr_ref_fai
    path chr_ref_dict

    output:
    path "happy_${caller}"

    script:
    """
    /opt/hap.py/bin/hap.py \\
      "${truth_vcf}" \\
      "${query_vcf}" \\
      -f "${truth_bed}" \\
      -r "${chr_ref}" \\
      -o "happy_${caller}" \\
      --engine=vcfeval \\
      --threads=${params.happy_threads}
    """
}

workflow {
    bam_file = file(params.bam)
    bai_file = file(params.bam + ".bai")
    ref_file = file(params.ref)
    ref_fai = file(params.ref + ".fai")
    ref_dict = file(params.ref.replace(".fna", ".dict").replace(".fasta", ".dict"))
    bed_file = file(params.chr_bed)
    truth_vcf = file(params.truth_vcf)
    truth_bed = file(params.truth_bed)

    // Generate mapping from accession to chromosome name
    mapping_file = file("accession_to_chr.map")
    file(mapping_file).text = ref_fai
        .readLines()
        .findAll { it.startsWith('NC_0000') }
        .collect { line ->
            def parts = line.split('\t')
            def acc = parts[0]
            def num = (acc =~ /NC_0000(\d+)\./)[0][1] as int
            if (num >= 1 && num <= 22) {
                return "${acc} chr${num}"
            } else {
                return null
            }
        }
        .findAll { it != null }
        .join('\n')

    // Prepare chr-prefixed reference for hap.py
    (chr_ref, chr_ref_fai, chr_ref_dict) = prepareChrRef(ref_file, mapping_file)

    // Run variant callers
    deep_vcf = deepvariant(bam_file, bai_file, ref_file, ref_fai, ref_dict, bed_file)
    clair_vcf = clair3(bam_file, bai_file, ref_file, ref_fai, bed_file)

    // Create a channel of raw VCFs with caller labels
    deep_ch = deep_vcf.map { vcf -> tuple('deepvariant', vcf) }
    clair_ch = clair_vcf.map { vcf -> tuple('clair3', vcf) }
    vcf_raw_ch = deep_ch.mix(clair_ch)

    // Rename contigs
    vcf_renamed_ch = renameContigs(vcf_raw_ch, mapping_file)

    // Benchmark using chr-prefixed reference
    hapPy(vcf_renamed_ch, truth_vcf, truth_bed, chr_ref, chr_ref_fai, chr_ref_dict)
}
