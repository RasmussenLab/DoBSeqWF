process MERGEBAM {
    label 'process_medium'
    tag "Merge bam files - $sample_id"
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*"), path(ubam_file, stageAs: "raw/*")
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_raw.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk MergeBamAlignment              \
        TMP_DIR=.                       \
        UNMAPPED=${ubam_file}           \
        ALIGNED=${bam_file}             \
        O=${sample_id}_raw.bam          \
        R=${db}                         \
        CLIP_ADAPTERS=false             \
        VALIDATION_STRINGENCY=SILENT    \
        CREATE_INDEX=true               \
        EXPECTED_ORIENTATIONS=FR        \
        MAX_GAPS=-1                     \
        SO=coordinate                   \
        ALIGNER_PROPER_PAIR_FLAGS=false
    """
    
    stub:
    """
    touch "${sample_id}_raw.bam"
    """
}