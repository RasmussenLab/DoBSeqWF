process FLAGSTAT {
    tag "Flagstat - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path("${sample_id}.flagstat"), emit: flagstat

    script:
    """
    samtools flagstat   \
        ${bam_file}     \
        > "${sample_id}.flagstat"
    """
    
    stub:
    """
    touch "${sample_id}.flagstat"
    """
}

