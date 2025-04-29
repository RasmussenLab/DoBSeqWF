process GROUP_UMI {
    label 'process_medium'
    tag "Group reads by UMI - $sample_id"

    conda "$projectDir/envs/fgbio/environment.yaml"
    container params.container.fgbio

    publishDir "${params.outputDir}/log/group_umi/", pattern: "${sample_id}.family_size_histogram.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file, stageAs: 'raw/*')

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: grouped_bam_file
    path "${sample_id}.family_size_histogram.txt", emit: family_metrics

    script:
    """
    ${params.fgbio} --tmp-dir=. --compression 1 --async-io GroupReadsByUmi           \
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