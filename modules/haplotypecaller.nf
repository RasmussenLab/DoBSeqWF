process HAPLOTYPECALLER {
    tag "HaplotypeCaller - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/haplotypecaller/", pattern: "${sample_id}.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/haplotypecaller/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.GATK.vcf.gz"), emit: gatk_vcf_file
    path "${sample_id}.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk HaplotypeCaller                                \
        --max-reads-per-alignment-start 0               \
        -R ${db}                                        \
        -I ${bam_file}                                  \
        -L ${bedfile}                                   \
        --disable-read-filter NotDuplicateReadFilter    \
        -O "${sample_id}.GATK.vcf.gz"                   \
        -ploidy ${params.ploidy}                        \
        --max-alternate-alleles 3                       \
        --max-num-haplotypes-in-population 1000         \
        > >(tee -a "${sample_id}.haplotypecaller.log")  \
        2> >(tee -a "${sample_id}.haplotypecaller.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.GATK.vcf.gz" "${sample_id}.haplotypecaller.log"
    """
}

