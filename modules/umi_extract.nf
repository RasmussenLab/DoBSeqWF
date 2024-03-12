process UMI_EXTRACT {
    tag "Extract UMI from uBAM - $sample_id"
    // Convert FastQ files to unaligned BAM files.
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*")

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: umi_extracted_ubam_file

    script:
    """
    fgbio ExtractUmisFromBam                \
        --input=${bam_file}                 \
        --output=${sample_id}.bam           \
        --read-structure=5M2S+T 5M2S+T      \
        --molecular-index-tags=ZA ZB        \
        --single-tag=RX
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}