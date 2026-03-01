#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * VARIANT CALLING PIPELINE
 * -------------------------
 * Input: FASTQ of PacBio HiFi reads (or pre‑aligned BAM with skip_alignment=true).
 * Steps:
 *   1. If needed, align reads with minimap2 (map-hifi) → sorted BAM + index.
 *   2. Prepare reference indexes (.fai with samtools, .dict with Picard).
 *   3. Generate chr1-22 BED from reference.
 *   4. Run DeepVariant (using run_deepvariant wrapper) and Clair3.
 *   5. Rename VCF contigs from accessions (NC_0000XX.X) to chr1..22.
 *   6. Create a chr‑prefixed reference for hap.py (two processes: rename/index, then dictionary).
 *   7. Benchmark each renamed VCF with hap.py (vcfeval engine).
 * Output: Benchmark summary CSVs in results/benchmark/.
 */

params.baseDir          = "/home/kashif/nazar/acb"
params.fastq            = "${params.baseDir}/m21009_241011_231051.hifi_reads.fastq.gz"
params.bam              = "${params.baseDir}/aligned_reads.bam"      // used only if skip_alignment=true
params.ref              = "${params.baseDir}/GCF_000001405.40_GRCh38.p14_genomic.fna"
params.truthVcf         = "${params.baseDir}/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truthBed         = "${params.baseDir}/giab_truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
params.outdir           = "${params.baseDir}/results"
params.skip_alignment   = false   // set true if aligned BAM already exists

// Container images
params.samtools_img     = "docker://staphb/samtools:latest"
params.minimap2_img     = "docker://staphb/minimap2:latest"
params.picard_img       = "docker://broadinstitute/picard:latest"
params.deepvariant_img  = "docker://google/deepvariant:1.6.0"
params.clair3_img       = "docker://hkubal/clair3:latest"
params.bcftools_img     = "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
params.happy_img        = "docker://quay.io/biocontainers/hap.py:0.3.7--py27_0"

// Resource settings
params.num_cpus         = 8

// -----------------------------------------------------------
// Processes
// -----------------------------------------------------------

process alignReads {
    tag "alignReads"
    container params.minimap2_img
    publishDir "${params.outdir}/bam", mode: 'copy'
    errorStrategy 'terminate'
    cpus params.num_cpus

    input:
    path fastq
    path ref
    path ref_fai

    output:
    tuple path("aligned_reads.bam"), path("aligned_reads.bam.bai")

    script:
    """
    minimap2 -ax map-hifi -t ${task.cpus} ${ref} ${fastq} \
        | samtools sort -@ ${task.cpus} -o aligned_reads.bam -
    samtools index aligned_reads.bam
    """
}

process indexFasta {
    tag "indexFasta"
    container params.samtools_img
    publishDir "${params.baseDir}", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path ref

    output:
    path "*.fai"

    script:
    """
    samtools faidx ${ref}
    """
}

process dictFasta {
    tag "dictFasta"
    container params.picard_img
    publishDir "${params.baseDir}", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path ref

    output:
    path "*.dict"

    script:
    def dictFile = ref.name.replace('.fna','.dict').replace('.fasta','.dict')
    """
    picard CreateSequenceDictionary R=${ref} O=${dictFile}
    """
}

process makeBed {
    tag "makeBed"
    container params.samtools_img
    publishDir "${params.baseDir}", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path ref_fai

    output:
    path "chr1-22.bed"

    script:
    """
    grep -E '^NC_0000(0[1-9]|1[0-9]|2[0-2])\\.' ${ref_fai} \
        | cut -f1,2 | awk -v OFS='\\t' '{print \$1, 0, \$2}' > chr1-22.bed
    """
}

process deepvariant {
    tag "deepvariant"
    container params.deepvariant_img
    publishDir "${params.outdir}/deepvariant", mode: 'copy'
    errorStrategy 'terminate'
    cpus params.num_cpus
    memory '16 GB'

    input:
    tuple path(bam), path(bai), path(ref), path(fai), path(dict), path(bed)

    output:
    path "deepvariant.vcf.gz"

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --regions ${bed} \
        --output_vcf deepvariant.vcf.gz \
        --num_shards ${task.cpus}
    """
}

process clair3 {
    tag "clair3"
    container params.clair3_img
    publishDir "${params.outdir}/clair3", mode: 'copy'
    errorStrategy 'terminate'
    cpus params.num_cpus
    memory '16 GB'

    input:
    tuple path(bam), path(bai), path(ref), path(fai), path(bed)

    output:
    path "clair3.vcf.gz"

    script:
    """
    mkdir -p clair3_tmp
    run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --threads=${task.cpus} \
        --platform="hifi" \
        --model_path="/opt/models/hifi" \
        --bed_fn=${bed} \
        --output=clair3_tmp
    cp clair3_tmp/merge_output.vcf.gz clair3.vcf.gz
    """
}

process makeMapping {
    tag "makeMapping"
    container params.samtools_img
    publishDir "${params.baseDir}", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path ref_fai

    output:
    path "accession_to_chr.map"

    script:
    """
    awk '\$1 ~ /^NC_0000/ {split(\$1,a,"."); n=a[1]; sub(/NC_0000/,"",n); n=int(n); if(n>=1 && n<=22) print \$1, "chr"n}' ${ref_fai} > accession_to_chr.map
    """
}

process renameContigs {
    tag "renameContigs-${caller}"
    container params.bcftools_img
    publishDir "${params.outdir}/renamed_vcfs", mode: 'copy'
    errorStrategy 'terminate'

    input:
    tuple val(caller), path(vcf)
    path mapping

    output:
    tuple val(caller), path("${caller}.renamed.vcf.gz")

    script:
    """
    bcftools annotate --rename-chrs ${mapping} -Oz -o ${caller}.renamed.vcf.gz ${vcf}
    bcftools index ${caller}.renamed.vcf.gz
    """
}

// Process to rename reference headers and index (samtools container)
process renameRef {
    tag "renameRef"
    container params.samtools_img
    publishDir "${params.outdir}/chr_reference", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path ref
    path mapping

    output:
    tuple path("chr_ref.fna"), path("chr_ref.fna.fai")

    script:
    """
    awk 'NR==FNR {map[\$1]=\$2; next}
         /^>/ {
             split(\$0,a," ");
             h=substr(a[1],2);
             if(h in map) {
                 print ">"map[h] " " substr(\$0, index(\$0,\$2));
             } else {
                 print \$0
             }
             next
         }
         { print }' ${mapping} ${ref} > chr_ref.fna
    samtools faidx chr_ref.fna
    """
}

// Process to create dictionary for renamed reference (Picard container)
process dictRenamedRef {
    tag "dictRenamedRef"
    container params.picard_img
    publishDir "${params.outdir}/chr_reference", mode: 'copy'
    errorStrategy 'terminate'

    input:
    path chr_ref

    output:
    path "chr_ref.dict"

    script:
    """
    picard CreateSequenceDictionary R=${chr_ref} O=chr_ref.dict
    """
}

process hapPy {
    tag "hapPy-${caller}"
    container params.happy_img
    publishDir "${params.outdir}/benchmark", mode: 'copy'
    errorStrategy 'terminate'
    cpus params.num_cpus
    memory '8 GB'

    input:
    tuple val(caller), path(query_vcf)
    path truth_vcf
    path truth_bed
    path chr_ref
    path chr_ref_fai
    path chr_ref_dict

    output:
    path "happy_${caller}*"

    script:
    """
    /usr/local/bin/hap.py \
        ${truth_vcf} \
        ${query_vcf} \
        -f ${truth_bed} \
        -r ${chr_ref} \
        -o happy_${caller} \
        --engine=vcfeval \
        --threads=${task.cpus}
    """
}

// -----------------------------------------------------------
// Workflow
// -----------------------------------------------------------
workflow {

    // --------------------------------------------------
    // Input channels
    // --------------------------------------------------
    ref_ch        = Channel.fromPath(params.ref)
    truth_vcf_ch  = Channel.fromPath(params.truthVcf)
    truth_bed_ch  = Channel.fromPath(params.truthBed)

    // --------------------------------------------------
    // Reference indexes
    // --------------------------------------------------
    fai_ch  = indexFasta(ref_ch)
    dict_ch = dictFasta(ref_ch)

    // Combine fai and dict into a single tuple for later use
    ref_indexes = fai_ch.combine(dict_ch)

    // BED from FAI
    bed_ch = makeBed(fai_ch)

    // --------------------------------------------------
    // Reads channel
    // --------------------------------------------------
    if (params.skip_alignment) {
        // Use provided BAM and its index
        reads_ch = Channel.fromPath(params.bam)
            .map { bam -> tuple(bam, file(bam + ".bai")) }
    } else {
        // Align from FASTQ
        fastq_ch = Channel.fromPath(params.fastq)
        reads_ch = alignReads(fastq_ch, ref_ch, fai_ch)
    }

    // --------------------------------------------------
    // Build DeepVariant input: (bam, bai, ref, fai, dict, bed)
    // --------------------------------------------------
    dv_input = reads_ch
        .combine(ref_ch)
        .combine(ref_indexes)
        .combine(bed_ch)
        .map { bam, bai, ref, fai, dict, bed ->
            tuple(bam, bai, ref, fai, dict, bed)
        }

    dv_vcf = deepvariant(dv_input)

    // --------------------------------------------------
    // Clair3 input: (bam, bai, ref, fai, bed) – drop dict
    // --------------------------------------------------
    cl_input = dv_input.map { bam, bai, ref, fai, dict, bed ->
        tuple(bam, bai, ref, fai, bed)
    }

    cl_vcf = clair3(cl_input)

    // --------------------------------------------------
    // Contig renaming
    // --------------------------------------------------
    mapping_ch = makeMapping(fai_ch)

    dv_label = dv_vcf.map { vcf -> tuple("deepvariant", vcf) }
    cl_label = cl_vcf.map { vcf -> tuple("clair3", vcf) }

    all_vcfs = dv_label.mix(cl_label)

    renamed_vcfs = renameContigs(all_vcfs, mapping_ch)

    // --------------------------------------------------
    // Prepare chr-prefixed reference for hap.py
    // Two steps: rename & index (samtools), then dictionary (Picard)
    // --------------------------------------------------
    // Step 1: rename reference and create .fai
    chr_ref_tup = renameRef(ref_ch, mapping_ch)                 // emits (chr_ref.fna, chr_ref.fna.fai)

    // Step 2: create dictionary from renamed reference
    chr_dict = dictRenamedRef( chr_ref_tup.map { ref, fai -> ref } )  // emits chr_ref.dict

    // Step 3: combine correctly (using combine, not join)
    chr_full = chr_ref_tup
        .combine(chr_dict)
        .map { ref, fai, dict ->
            tuple(ref, fai, dict)
        }

    // Split into separate channels for downstream use
    chr_ref      = chr_full.map { ref, fai, dict -> ref }
    chr_ref_fai  = chr_full.map { ref, fai, dict -> fai }
    chr_ref_dict = chr_full.map { ref, fai, dict -> dict }

    // --------------------------------------------------
    // hap.py input: combine renamed_vcfs with truth and chr reference
    // --------------------------------------------------
    hap_input = renamed_vcfs
        .combine(truth_vcf_ch)
        .combine(truth_bed_ch)
        .combine(chr_ref)
        .combine(chr_ref_fai)
        .combine(chr_ref_dict)
        .map { caller, query_vcf, truth_vcf, truth_bed, ref, fai, dict ->
            tuple(caller, query_vcf, truth_vcf, truth_bed, ref, fai, dict)
        }

    hapPy(hap_input)
}