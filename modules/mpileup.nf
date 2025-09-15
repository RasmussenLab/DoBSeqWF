process MPILEUP {
    label 'process_low'
    tag "BAM->mpileup - $sample_id"
    
    conda "$projectDir/envs/samtools/environment.yaml"
    container params.container.samtools

    publishDir "${params.outputDir}/mpileup/", mode:'copy', pattern: "${sample_id}.mpileup.gz"

    input:
    tuple val(sample_id), path(bam_file), path(index)
    path(reference)
    path(bedfile)

    output:
    path("${sample_id}.mpileup.gz"), emit: mpileup_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    samtools mpileup \
        -l ${bedfile} \
        --max-depth 50000 \
        --min-BQ 0 \
        --output-MQ \
        --fasta-ref ${db} \
        ${bam_file} \
            | gzip > ${sample_id}.mpileup.gz
    """
    
    stub:
    """
    touch "${sample_id}.mpileup.gz"
    """
}

