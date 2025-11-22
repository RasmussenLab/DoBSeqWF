process BAM {
    label 'process_low'
    tag "CRAM->BAM - $sample_id"
    
    conda "$projectDir/envs/samtools/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.samtools : params.container.docker.samtools

    input:
    tuple val(sample_id), path(cram_file)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_file

    script:
    def db = file(params.reference_genome).getName() + ".fasta"
    """
    samtools view                       \
        -@ ${task.cpus-1}               \
        -T ${db}                        \
        -b                              \
        -o "${sample_id}.bam"          \
        ${cram_file}

    ## -b: output in bam format
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    """
}

