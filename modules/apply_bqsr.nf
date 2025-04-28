process APPLY_BQSR {
    tag "Apply score recalibration - $sample_id"

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*"), path(bqsr_table)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: corrected_bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk ApplyBQSR              \
        -R ${db}                \
        -I ${bam_file}          \
        -bqsr ${bqsr_table}     \
        -O "${sample_id}.bam"
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}
