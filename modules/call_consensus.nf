process CALL_CONSENSUS {
    tag "Call consensus reads - $sample_id"

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.ubam"), emit: consensus_ubam_file

    script:
    """
    fgbio --compression 1 --async-io CallDuplexConsensusReads          \
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