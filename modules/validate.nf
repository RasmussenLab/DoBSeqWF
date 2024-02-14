process VALIDATE {
    tag "Validate bam file - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/validation/", pattern: "${sample_id}.validation.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.validation.log"), emit: validation_log

    script:
    """
    gatk ValidateSamFile                \
        -I ${bam_file}                  \
        -M SUMMARY                      \
        --IGNORE_WARNINGS true          \
        -O "${sample_id}.validation.log"
    """
    
    stub:
    """
    touch "${sample_id}.validation.log"
    """
}

