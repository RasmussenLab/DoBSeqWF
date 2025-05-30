process CALL_CONSENSUS {
    label 'process_low'
    tag "Call consensus reads - $sample_id"

    conda "$projectDir/envs/fgbio/environment.yaml"
    container params.container.fgbio

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.ubam"), emit: consensus_ubam_file

    script:
    """
    ${params.fgbio} --tmp-dir=. --compression 1 --async-io CallDuplexConsensusReads          \
        --input=${bam_file}                 \
        --output="${sample_id}.ubam"        \
        --error-rate-pre-umi=45             \
        --error-rate-post-umi=30            \
        --min-input-base-quality=30         \
        --min-reads 1 0 0
    """
    
    stub:
    """
    touch "${sample_id}.ubam"
    """
}