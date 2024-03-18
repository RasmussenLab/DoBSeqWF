process GROUP_UMI {
    tag "Group reads by UMI - $sample_id"

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/group_umi/", pattern: "${sample_id}.family_size_histogram.txt", mode:'copy'


    input:
    tuple val(sample_id), path(bam_file, stageAs: 'raw/*')

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: grouped_bam_file
    path "${sample_id}.family_size_histogram.txt"

    script:
    """
    fgbio --compression 1 --async-io GroupReadsByUmi           \
        --strategy=paired           \
        --input=${bam_file}         \
        --output="${sample_id}.bam" \
        --raw-tag=RX                \
        --min-map-q=10              \
        --edits=1                   \
        --family-size-histogram="${sample_id}.family_size_histogram.txt"
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    touch "${sample_id}.family_size_histogram.txt"
    """
}