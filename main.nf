#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// —— PARAMETERS ——
params.data            = params.data            ?: '/Data/file'
params.outdir          = params.outdir          ?: '/results'
params.reference       = params.reference       ?: 'reference_H37RV.fasta'
params.primers         = params.primers         ?: 'tb-amplicon-primers.bed'
params.compare_script  = params.compare_script  ?: 'compare_mutations.py'

workflow {
    // —— CHANNEL SETUP ——
    primers_txt    = Channel.fromPath(params.primers).first()
    reference      = Channel.fromPath(params.reference).first()
    compare_script = Channel.fromPath(params.compare_script).first()
    fastq_ch       = Channel.fromPath("${params.data}/*.fastq.gz")
                         .ifEmpty { error "No FASTQ files found in: ${params.data}" }
    sample_ch      = fastq_ch.map { it.baseName.replaceFirst(/\.fastq(?:\.gz)?$/, '') }

    // —— PIPELINE STEPS ——
    qual_ch         = runQualityTrimming(fastq_ch)
    align_ch        = runAlignment(qual_ch, reference)
    sorted_ch       = runSortAndIndex(align_ch)
    variant_ch      = runVariantCalling(sorted_ch, reference)
    filtered_vcf_ch = runFilterVariants(variant_ch)
    mutations_ch    = runConvertToTSV(filtered_vcf_ch)

    compareMutations(
        mutations_ch,
        sample_ch,
        Channel.value(params.outdir),
        reference,
        compare_script
    )
}

// —— PROCESSES ——

process runQualityTrimming {
    input:
      path fastq
    output:
      path "quality_trimmed_${fastq.baseName}.fastq.gz"
    script:
    """
    seqkit rename ${fastq} -o renamed.fastq
    filtlong --min_length 10 --keep_percent 90 renamed.fastq | gzip > quality_trimmed_${fastq.baseName}.fastq.gz
    """
}

process runAlignment {
    input:
      path trimmed
      path reference
    output:
      path "aligned_${trimmed.baseName}.sam"
    script:
    """
    bwa index -p goober ${reference}
    bwa mem -P goober ${trimmed} > aligned_${trimmed.baseName}.sam
    """
}

process runSortAndIndex {
    input:
      path sam
    output:
      path "aligned_sorted_${sam.baseName}.bam"
    script:
    """
    samtools sort ${sam} -o aligned_sorted_${sam.baseName}.bam
    samtools index aligned_sorted_${sam.baseName}.bam
    """
}

process runTrimmingIvar {
    input:
      path bam
      path primers_txt
    output:
      path "trimmed_${bam.baseName}.bam"
    script:
    """
    ivar trim -b ${primers_txt} -i ${bam} -p trimmed_${bam.baseName}.bam -q 2 -x 1000
    """
}

process runVariantCalling {
    input:
      path bam
      path reference
    output:
      path "variants_${bam.baseName}.vcf"
    script:
    """
    freebayes -p 1 -f ${reference} ${bam} > variants_${bam.baseName}.vcf
    """
}

process runFilterVariants {
    input:
      path vcf
    output:
      path "filtered_${vcf.baseName}.vcf"
    script:
    """
    bcftools view -i \\
      'FMT/GT="1" && QUAL>=20 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0.9' \\
      ${vcf} > filtered_${vcf.baseName}.vcf
    """
}

process runConvertToTSV {
    input:
      path vcf
    output:
      path "variants_${vcf.baseName}.tsv"
    script:
    """
    (
      echo -e "CHROM\tPOS\tALT\tREF\tQUAL\tINFO"
      bcftools query -f '%CHROM\\t%POS\\t%ALT\\t%REF\\t%QUAL\\t%INFO\\n' ${vcf}
    ) > variants_${vcf.baseName}.tsv
    """
}

process compareMutations {
    input:
      path mutations
      val sample
      val outdir
      path reference
      path script
    script:
    """
    python3 ${script} ${mutations} ${sample} ${outdir} ${reference}
    """
}




