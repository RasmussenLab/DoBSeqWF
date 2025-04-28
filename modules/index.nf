process INDEX {
    tag "Index bam file - $sample_id"

    conda "$projectDir/envs/samtools/environment.yaml"
    container params.container.samtools

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path(bam_file), path("${bam_file}.bai"), emit: bam_file_w_index

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

