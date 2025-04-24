process EXTRACT_UMI {
    tag "Extract UMI from uBAM - $sample_id"
    // Convert FastQ files to unaligned BAM files.
    
    // UMI bases are stored in three tags: 
    //  ZA for the 5' UMI
    //  ZB for the 3' UMI
    //  RX for a string containing both the 5' and 3' UMIs separated by a dash.
    // Twist UMI adapters are 5 base pairs with a 2 base pair skip, resulting in the read structure 5M2S+T

    conda "$projectDir/envs/fgbio/environment.yaml"

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*")

    output:
    tuple val(sample_id), path("${sample_id}.ubam"), emit: umi_extracted_ubam_file

    script:
    """
    ${params.fgbio} --tmp-dir=. --compression 1 --async-io ExtractUmisFromBam                \
        --input=${bam_file}                 \
        --output=${sample_id}.ubam           \
        --read-structure=5M2S+T 5M2S+T      \
        --molecular-index-tags=ZA ZB        \
        --single-tag=RX
    """
    
    stub:
    """
    touch "${sample_id}.ubam"
    """
}