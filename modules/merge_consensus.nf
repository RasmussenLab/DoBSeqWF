process MERGE_CONSENSUS {
    label 'process_medium'
    tag "Merge consensus bam files - $sample_id"
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*"), path(ubam_file, stageAs: "raw/*")
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk --java-options -Xmx20g MergeBamAlignment              \
        TMP_DIR=.                       \
        UNMAPPED=${ubam_file}           \
        ALIGNED=${bam_file}             \
        O=${sample_id}.bam              \
        R=${db}                         \
        CLIP_ADAPTERS=false             \
        VALIDATION_STRINGENCY=SILENT    \
        CREATE_INDEX=true               \
        EXPECTED_ORIENTATIONS=FR        \
        MAX_GAPS=-1                     \
        SO=coordinate                   \
        ALIGNER_PROPER_PAIR_FLAGS=false \
        ATTRIBUTES_TO_RETAIN=X0         \
        ATTRIBUTES_TO_RETAIN=ZS         \
        ATTRIBUTES_TO_RETAIN=ZI         \
        ATTRIBUTES_TO_RETAIN=ZM         \
        ATTRIBUTES_TO_RETAIN=ZC         \
        ATTRIBUTES_TO_RETAIN=ZN         \
        ATTRIBUTES_TO_RETAIN=ad         \
        ATTRIBUTES_TO_RETAIN=bd         \
        ATTRIBUTES_TO_RETAIN=cd         \
        ATTRIBUTES_TO_RETAIN=ae         \
        ATTRIBUTES_TO_RETAIN=be         \
        ATTRIBUTES_TO_RETAIN=ce
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}