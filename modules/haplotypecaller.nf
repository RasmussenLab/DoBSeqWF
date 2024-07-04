process HAPLOTYPECALLER {
    tag "HaplotypeCaller - $sample_id - $interval"
    // Call variants using GATK - HaplotypeCaller
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }
    publishDir "${params.outputDir}/log/haplotypecaller/${sample_id}/", pattern: "${sample_id}*.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file), val(interval)
    path reference_genome
    path interval_list

    output:
    tuple val(sample_id), path("${sample_id}*.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}_vcf.tsv", emit: vcf_info
    path "${sample_id}*.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def id = interval ? "${sample_id}.${interval}" : sample_id
    def add_intervals = interval ? "-L ${interval}" : ""
    def add_intersection = interval ? "--interval-set-rule INTERSECTION" : ""
    def publishDir = file(params.outputDir + "/variants/" + sample_id + ".GATK.vcf.gz")
    """
    gatk HaplotypeCaller                                \
        --max-reads-per-alignment-start 0               \
        -R ${db}                                        \
        -I ${bam_file}                                  \
        -L ${interval_list}                             \
        ${add_intervals}                                \
        ${add_intersection}                             \
        -O "${id}.GATK.vcf.gz"       \
        -ploidy ${params.ploidy}                        \
        --max-alternate-alleles 3                       \
        > >(tee -a "${id}.haplotypecaller.log")  \
        2> >(tee -a "${id}.haplotypecaller.log" >&2)
    
    echo -e "${sample_id}\t${publishDir}\tGATK" > "${sample_id}_vcf.tsv"
    """

    stub:
    """
    touch "${id}.GATK.vcf.gz" "${id}.haplotypecaller.log"
    touch "${sample_id}_vcf.tsv"
    """
}
