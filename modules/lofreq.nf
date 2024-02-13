process LOFREQ {
    tag "Lofreq call - $sample_id"
    // Call variants using Lofreq - call parallel
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/lofreq/", pattern: "${sample_id}.lofreq.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.lofreq.vcf.gz"), emit: lofreq_vcf_file
    path "${sample_id}.lofreq.log"

    script:
    """
    lofreq call-parallel                    \
        --call-indels                       \
        --pp-threads ${task.cpus}           \
        --force-overwrite                   \
        --no-default-filter                 \
        --sig 1                             \
        --bonf 1                            \
        -f ${reference_genome}              \
        -l ${bedfile}                       \
        -o "${sample_id}.lofreq.vcf.gz"     \
        ${bam_file}                         \
        &> "${sample_id}.lofreq.log"
    """
    
    stub:
    """
    touch "${sample_id}.lofreq.vcf.gz" "${sample_id}.lofreq.log"
    """
}

