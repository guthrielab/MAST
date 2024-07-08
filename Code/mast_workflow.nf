#!/usr/bin/env nextflow
nextflow.enable.dsl=2

"""
pipeline input paramaters
"""

params.data = '/some/data/file'
params.output_dir = './results'
params.reference = './Data/GCF_000195955.2_ASM19595v2_genomic.fasta'
params.primers = './Data/primers.txt'


log.info"""\
        M T B   A M P L I C O N   T O O L
        =================================
        data: $params.data
        output_dir: $params.output_dir
        reference: $params.reference
        primers: $params.primers
        """

process runTrimming {
    input:
    file fastq
    file primers

    output:
    file 'trimmed_file.fastq.gz'

    script:
    """
    cutadapt_options=""
    options=(" -g" " -a")
    index=0
    while IFS= read -r line; do
         cutadapt_options+="\${options[index]} \$line"
         index=\$(( (index + 1) % \${#options[@]} ))
    done < "$primers"

    echo "Constructed cutadapt options: \$cutadapt_options" > test.txt

    cutadapt \$cutadapt_options -o trimmed_file.fastq.gz $fastq
    """
}

process runQualityTrimming {
    input:
    file fastq_trimmed

    output:
    file 'quality_trimmed.fastq.gz'

    script:
    """
    filtlong --min_length 10 --max_length 1000000 --keep_percent 90 --target_bases 400000000 $fastq_trimmed | gzip > 'quality_trimmed.fastq.gz'
    """

}
 
process runAllignment {
    input:
    file trimmed_file
    file reference

    output:
    file 'alligned.sam'

    script:
    """
    bwa index -p goober $reference
    bwa mem -P goober $trimmed_file > alligned.sam
    """
}

process runSortAndIndex {

    input:
    file alligned

    output:
    file 'alligned.sorted.bam'

    script:
    """
    samtools sort $alligned > alligned.sorted.bam
    samtools index alligned.sorted.bam
    """
}

process runVariantCalling {

    input:
    file alligned
    file reference
    

    output:
    file 'variants.vcf'

    script:
    """
    freebayes -f $reference $alligned > variants.vcf
    """
}

process runFilterVariants {

    input:
    file rawvariants
    
    output:
    file 'filteredvariants.vcf'
    
    script:
    """
    bcftools view --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' $rawvariants > filteredvariants.vcf
    """
}

process runConvertToTSV {

    input:
    file variants

    output:
    file "variants_h37ra_2.tsv"

    script:
    """
    (
    echo -e "CHROM\tPOS\tALT\tREF\tQUAL\tINFO" &&
    bcftools query -f '%CHROM\t%POS\t%ALT\t%REF\t%QUAL\t%INFO\n' $variants
    ) | awk -F'\t' 'BEGIN{OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6}' > variants_h37ra_2.tsv  
    """
}

process compareMutations {

    input:
    file mutations
    file data
    file output_dir
    file reference

    script:
    """
    python3 ./Code/compare.py ${mutations} ${data.simpleName} ${output_dir} ${reference}
    """
}

workflow {
    primers_txt = Channel.fromPath(params.primers)
    reference = Channel.fromPath(params.reference)

    qual_ch = runQualityTrimming(Channel.fromPath(params.data))


    allignment_ch = runAllignment(qual_ch, reference)
    sorted_ch = runSortAndIndex(allignment_ch)
    variant_ch = runVariantCalling(sorted_ch, reference)
    raw_variant_ch = runFilterVariants(variant_ch)
    mutations_ch = runConvertToTSV(variant_ch)
    compareMutations(mutations_ch, Channel.fromPath(params.data), Channel.fromPath(params.output_dir), Channel.fromPath(params.reference))
}
