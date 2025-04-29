process ADDREADGROUP {
    label 'process_single'
    tag "Add read group - $sample_id"

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*")

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    """
    gatk AddOrReplaceReadGroups     \
        --INPUT ${bam_file}   		\
        --OUTPUT "${sample_id}.bam" \
        --RGID 4                    \
        --RGLB lib1                 \
        --RGPL illumina             \
        --RGPU unit1                \
        --RGSM ${sample_id}
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}
