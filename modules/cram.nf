process CRAM {
    tag "BAM->CRAM - $sample_id"
    
    conda "$projectDir/envs/samtools/environment.yaml"

    publishDir "${params.outputDir}/cram/", pattern: "${sample_id}.cram", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: cram_file
    path("${sample_id}_cram.tsv"), emit: cram_info

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def publishDir = file(params.outputDir + "/cram/" + sample_id + ".cram")
    """
    samtools view                       \
        -@ ${task.cpus-1}               \
        -T ${db}                        \
        -C                              \
        -o "${sample_id}.cram"          \
        ${bam_file}
    
    echo -e "${sample_id}\t${publishDir}" > "${sample_id}_cram.tsv"
    ## -c: output in cram format
    """
    
    stub:
    """
    touch "${sample_id}.cram"
    touch "${sample_id}_cram.tsv"
    """
}

