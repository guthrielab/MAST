#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if (params.help) {
  helpMessage()
  exit 0
}

def helpMessage() {
  log.info """
  ====================================================
  MAST - Mycobacteria Amplicon Sequencing Tool
  ====================================================

  Usage:
  nextflow run MAST --input <input_dir> --outdir <output_dir> -profile <conda>[options]

  Required parameters:
  --input                       Input directory containing FASTQ files
  --outdir                      Output directory (default: results)
  -profile                      Signify what container to use (conda, docker, or singularity)
  --help                        Show this help message

  Optional parameters:
  --reference                   Reference genome for alignment (default: Mtb H37Rv)
  --primers                     BED file to trim primers (default: Amplicon primers MAST/tb-amplicon-primers.bed)
  --mutations_csv               Resistance mutations CSV (default: MAST/Data/all_resistant_variants.csv)
  --lineage_csv                 Lineage CSV (default: MAST/Data/Lineage.csv)
  --template_docx               Report template (default: MAST/Data/tNGS_Report_Template.docx)
  --patient_info_csv            Patient info CSV (default: MAST/Data/patient_info.csv)

  For more information, visit: https://github.com/guthrielab/MAST
  """
}

// —— PARAMETERS ——
params.input            = params.input            ?: 'MAST/Data/file'
params.outdir           = params.outdir           ?: 'MAST/results'
params.reference        = 'MAST/reference_H37RV.fasta'
params.primers          = 'MAST/tb-amplicon-primers.bed'
params.compare_script   = 'MAST/compare_mutations.py'
params.mutations_csv    = 'MAST/Data/all_resistant_variants.csv'
params.lineage_csv      = 'MAST/Data/Lineage.csv'
params.template_docx    = 'MAST/Data/tNGS_Report_Template.docx'
params.patient_info_csv = 'MAST/Data/patient_info.csv'
params.regions_bed      = 'MAST/regions.bed'

workflow {
    // —— CHANNEL SETUP ——
    primers_txt      = Channel.fromPath(params.primers).first()
    reference        = Channel.fromPath(params.reference).first()
    compare_script   = Channel.fromPath(params.compare_script).first()
    mutations_csv    = Channel.fromPath(params.mutations_csv).first()
    lineage_csv      = Channel.fromPath(params.lineage_csv).first()
    template_docx    = Channel.fromPath(params.template_docx).first()
    patient_info_csv = Channel.fromPath(params.patient_info_csv).first()
    regions_bed      = Channel.fromPath(params.regions_bed).first()

    fastq_ch = Channel
      .fromPath("${params.input}/*.fastq.gz")
      .ifEmpty { error "No FASTQ files found in: ${params.input}" }

    reads = fastq_ch.map { f -> tuple( f.baseName.replaceFirst(/\.fastq(?:\.gz)?$/, ''), f ) }

    // —— PIPELINE STEPS ——
    qual_ch         = runQualityTrimming(reads)
    align_ch        = runAlignment(qual_ch, reference)

    // sorted_ch emits: tuple val(id), path(bam), path(bai)
    // The BAI is carried forward so compareMutations can open the BAM
    // with pysam without needing to re-index inside the container.
    sorted_ch       = runSortAndIndex(align_ch)

    primertrim_ch     = runPrimerTrimming(sorted_ch, primers_txt)
    variant_ch      = runVariantCalling(primertrim_ch, reference)
    filtered_vcf_ch = runFilterVariants(variant_ch)
    mutations_ch    = runConvertToTSV(filtered_vcf_ch)

    // Join mutations TSV with sorted BAM+BAI on sample id.
    // Result: tuple val(id), path(tsv), path(bam), path(bai)
    mutations_with_bam = mutations_ch.join(sorted_ch)

    compareMutations(
        mutations_with_bam,
        reference,
        compare_script,
        mutations_csv,
        lineage_csv,
        template_docx,
        patient_info_csv,
        regions_bed
    )
}

// —— PROCESSES ——

process runQualityTrimming {
    container 'quay.io/idolawoye/mast_quality-trimming'

    input:
      tuple val(id), path(fastq)
    output:
      tuple val(id), path('quality_trimmed_*.fastq.gz')
    script:
    """
    set -euo pipefail
    seqkit rename ${fastq} -o renamed.fastq
    filtlong --min_length 10 --keep_percent 90 renamed.fastq | gzip > quality_trimmed_${id}.fastq.gz
    """
}

process runAlignment {
    container 'quay.io/biocontainers/minimap2:2.31--h118bc1c_0'

    input:
      tuple val(id), path(trimmed)
      path(reference)
    output:
      tuple val(id), path('aligned_*.sam')
    script:
    """
    set -euo pipefail

    minimap2 -t ${task.cpus} -a ${reference} quality_trimmed_${id}.fastq.gz > aligned_${id}.sam

    """
}

process runSortAndIndex {
    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
      tuple val(id), path(sam)
    output:
      // Emit both BAM and its index so downstream processes can stage
      // them together — pysam needs the .bai next to the .bam
      tuple val(id), path('aligned_sorted_*.bam'), path('aligned_sorted_*.bam.bai')
    script:
    """
    set -euo pipefail
    samtools sort -O bam -m 2G ${sam} -o aligned_sorted_${id}.bam
    samtools index aligned_sorted_${id}.bam
    """
}

process runPrimerTrimming {
    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
      tuple val(id), path(bam), path(bai)
      path(primers_txt)
    output:
      tuple val(id), path('trimmed_*.bam')
    script:
    """
    set -euo pipefail
    samtools ampliconclip -b ${primers_txt} ${bam} -o trimmed_${id}.bam
    """
}

process runVariantCalling {
    container 'quay.io/biocontainers/freebayes:1.3.9--hbefcdb2_1'

    input:
      tuple val(id), path(bam)
      path(reference)
    output:
      tuple val(id), path('variants_*.vcf')
    script:
    """
    set -euo pipefail
    freebayes -p 1 -f ${reference} ${bam} > variants_${id}.vcf
    """
}

process runFilterVariants {
    container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'

    input:
      tuple val(id), path(vcf)
    output:
      tuple val(id), path('filtered_*.vcf')
    script:
    """
    set -euo pipefail
    bcftools view -i \
      'FMT/GT="1" && QUAL>=20 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.9' \
      ${vcf} > filtered_${id}.vcf
    """
}

process runConvertToTSV {
    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    input:
      tuple val(id), path(vcf)
    output:
      tuple val(id), path('variants_*.tsv')
    script:
    """
    set -euo pipefail
    (
      echo -e "CHROM\\tPOS\\tALT\\tREF\\tQUAL\\tINFO"
      bcftools query -f '%CHROM\\t%POS\\t%ALT\\t%REF\\t%QUAL\\t%INFO\\n' ${vcf}
    ) > variants_${id}.tsv
    """
}

process compareMutations {
    container 'quay.io/idolawoye/mast_compare-mutations:v1.1.1'

    publishDir "${params.outdir}", mode: 'copy'

    input:
      // BAM and BAI are staged into the work directory by Nextflow.
      // pysam will find the .bai automatically next to the .bam.
      tuple val(id), path(mutations), path(bam), path(bai)
      path(reference)
      path(script)
      path(mutations_csv)
      path(lineage_csv)
      path(template_docx)
      path(patient_info_csv)
      path(regions_bed)
    output:
      path("${id}_report.docx")
      path("${id}_results.tsv")
    script:
    """
    set -euo pipefail
    python3 ${script} \
        ${mutations} \
        ${id} \
        . \
        ${reference} \
        ${mutations_csv} \
        ${lineage_csv} \
        ${template_docx} \
        ${patient_info_csv} \
        ${bam} \
        ${regions_bed}
    """
}
