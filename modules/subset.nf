process SUBSET {
    tag "Subset alignment and convert to bam - $sample_id"

    conda "$projectDir/envs/samtools/environment.yaml"

    input:
    tuple val(sample_id), path(bam_file, stageAs: 'raw/*')
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    samtools view               \
        -@ ${task.cpus}         \
        -b                      \
        -L ${bedfile}           \
        -o ${sample_id}.bam     \
        -T ${db}                \
        ${bam_file}
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}

