process LOFREQ {
    tag "Lofreq call - $sample_id"
    // Call variants using Lofreq - call parallel
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/lofreq/", pattern: "${sample_id}.lofreq.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.lofreq.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.lofreq.vcf.gz"), emit: vcf_file
    path "${sample_id}_vcf.tsv", emit: vcf_info
    path "${sample_id}.lofreq.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def publishDir = file(params.outputDir + "/variants/" + sample_id + ".lofreq.vcf.gz")
    """
    lofreq call-parallel                            \
        --call-indels                               \
        --pp-threads ${task.cpus}                   \
        --force-overwrite                           \
        --no-default-filter                         \
        --sig 1                                     \
        --bonf 1                                    \
        -f ${db}                                    \
        -l ${bedfile}                               \
        -o "${sample_id}.lofreq.vcf.gz"             \
        ${bam_file}                                 \
        > >(tee -a "${sample_id}.lofreq.log")       \
        2> >(tee -a "${sample_id}.lofreq.log" >&2)

    echo -e "${sample_id}\t${publishDir}\tlofreq" > "${sample_id}_vcf.tsv"
    """
    
    stub:
    """
    touch "${sample_id}.lofreq.vcf.gz"
    touch "${sample_id}.lofreq.log"
    touch "${sample_id}_vcf.tsv"
    """
}

