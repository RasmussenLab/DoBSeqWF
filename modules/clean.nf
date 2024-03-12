process CLEAN {
    tag "Clean bam file - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.clean.bam"), emit: clean_bam_file

    script:
    """
    samtools view                       \
        -F 1024                         \
        -F 512                          \
        -F 2048                         \
        -q 20                           \
        -b                              \
        -o "${sample_id}.clean.bam"     \
        ${bam_file}

    ## -F 1024: exclude read that fails platform/vendor quality checks
    ## -F 512: exclude read that is PCR or optical duplicate
    ## -F 2048: exclude supplementary alignments
    ## -q 20: exclude reads with mapping quality less than 20
    ## -b: output in BAM format
    """
    
    stub:
    """
    touch "${sample_id}.clean.bam"
    """
}

