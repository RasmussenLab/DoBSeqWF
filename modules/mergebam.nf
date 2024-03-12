process MERGEBAM {
    tag "Merge bam files - $sample_id"
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file, stageAs: "raw/*"), path(ubam_file, stageAs: "raw/*")
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk MergeBamAlignment              \
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
        ALIGNER_PROPER_PAIR_FLAGS=false
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}