process INDEX {
    tag "Index bam file - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${bam_file}.bai"), emit: bam_index_file

    script:
    """
    samtools index   \
        ${bam_file}
    """
    
    stub:
    """
    touch "${bam_file}.bai"
    """
}

