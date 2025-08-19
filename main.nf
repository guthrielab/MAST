#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// —— PARAMETERS ——
params.data            = params.data            ?: 'MAST/Data/file'
params.outdir          = params.outdir          ?: 'MAST/results'
params.reference       = params.reference       ?: 'MAST/reference_H37RV.fasta'
params.primers         = params.primers         ?: 'MAST/tb-amplicon-primers.bed'
params.compare_script  = params.compare_script  ?: 'MAST/compare_mutations.py'

workflow {
    // —— CHANNEL SETUP ——
    primers_txt    = Channel.fromPath(params.primers).first()
    reference      = Channel.fromPath(params.reference).first()
    compare_script = Channel.fromPath(params.compare_script).first()

    fastq_ch = Channel
      .fromPath("${params.data}/*.fastq.gz")
      .ifEmpty { error "No FASTQ files found in: ${params.data}" }


    reads = fastq_ch.map { f -> tuple( f.baseName.replaceFirst(/\.fastq(?:\.gz)?$/, ''), f ) }

    // —— PIPELINE STEPS ——
    qual_ch         = runQualityTrimming(reads)
    align_ch        = runAlignment(qual_ch, reference)
    sorted_ch       = runSortAndIndex(align_ch)
    // runTrimmingIvar exists but remains unused to preserve original behavior
    variant_ch      = runVariantCalling(sorted_ch, reference)
    filtered_vcf_ch = runFilterVariants(variant_ch)
    mutations_ch    = runConvertToTSV(filtered_vcf_ch)

    compareMutations(
        mutations_ch,
        Channel.value(params.outdir),
        reference,
        compare_script
    )
}

// —— PROCESSES ——

process runQualityTrimming {
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
    input:
      tuple val(id), path(trimmed)
      path(reference)
    output:
      tuple val(id), path('aligned_*.sam')
    script:
    """
    set -euo pipefail
    bwa index -p goober ${reference}
    bwa mem -P goober ${trimmed} > aligned_${id}.sam
    """
}

process runSortAndIndex {
    input:
      tuple val(id), path(sam)
    output:
      tuple val(id), path('aligned_sorted_*.bam')
    script:
    """
    set -euo pipefail
    samtools sort ${sam} -o aligned_sorted_${id}.bam
    samtools index aligned_sorted_${id}.bam
    """
}

process runTrimmingIvar {
    input:
      tuple val(id), path(bam)
      path(primers_txt)
    output:
      tuple val(id), path('trimmed_*.bam')
    script:
    """
    set -euo pipefail
    ivar trim -b ${primers_txt} -i ${bam} -p trimmed_${id}.bam -q 2 -x 1000
    """
}

process runVariantCalling {
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
    input:
      tuple val(id), path(mutations)
      val(outdir)
      path(reference)
      path(script)
    output:
      path('*_report.docx'), optional: true
    script:
    """
    set -euo pipefail
    mkdir -p ${outdir}
    python3 ${script} ${mutations} ${id} ${outdir} ${reference}
    """
}





