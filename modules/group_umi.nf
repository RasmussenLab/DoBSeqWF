process GROUP_UMI {
    tag "Group reads by UMI - $sample_id"

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file, stageAs: 'raw/*')

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: grouped_bam_file

    script:
    """
    fgbio GroupReadsByUmi           \
        --strategy=paired           \
        --input=${bam_file}         \
        --output="${sample_id}.bam" \
        --raw-tag=RX                \
        --min-map-q=10              \
        --edits=1
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}