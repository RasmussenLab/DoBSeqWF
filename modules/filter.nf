process FILTER {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/variants/filtered/", pattern: "${sample_id}.lofreq.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "input/*")
    val(alt_support)

    output:
    tuple val(sample_id), path("${sample_id}.lofreq.vcf.gz"), emit: filtered_vcf_file

    script:
    """
    bcftools view \
        --include 'INFO/DP4[2] + INFO/DP4[3] > ${alt_support}'  \
        --output "${sample_id}.lofreq.vcf.gz"                   \
        ${vcf_file}
    """

    stub:
    """
    touch "${sample_id}.lofreq.vcf.gz"
    """
}