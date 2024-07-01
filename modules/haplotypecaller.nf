process HAPLOTYPECALLER {
    tag "HaplotypeCaller - $sample_id - $interval"
    // Call variants using GATK - HaplotypeCaller
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }
    publishDir "${params.outputDir}/log/haplotypecaller/${sample_id}/", pattern: "${sample_id}.${interval}.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/intervals/", pattern: "${sample_id}.${interval}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file), val(interval)
    path reference_genome
    path interval_list

    output:
    tuple val(sample_id), path("${sample_id}*.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}*.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def add_intervals = interval ? "-L ${interval}" : ""
    def add_intersection = interval ? "--interval-set-rule INTERSECTION" : ""
    """
    gatk HaplotypeCaller                                \
        --max-reads-per-alignment-start 0               \
        -R ${db}                                        \
        -I ${bam_file}                                  \
        -L ${interval_list}                             \
        ${add_intervals}                                \
        ${add_intersection}                             \
        -O "${sample_id}.${interval}.GATK.vcf.gz"       \
        -ploidy ${params.ploidy}                        \
        --max-alternate-alleles 3                       \
        > >(tee -a "${sample_id}.${interval}.haplotypecaller.log")  \
        2> >(tee -a "${sample_id}.${interval}.haplotypecaller.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${interval}.GATK.vcf.gz" "${sample_id}.${interval}.haplotypecaller.log"
    """
}
