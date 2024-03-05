process BAM {
    tag "CRAM->BAM - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(cram_file)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    samtools view                       \
        -@ ${task.cpus-1}               \
        -T ${db}                        \
        -b                              \
        -o "${sample_id}.bam"          \
        ${cram_file}

    ## -b: output in bam format
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}

