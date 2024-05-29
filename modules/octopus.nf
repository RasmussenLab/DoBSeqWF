process OCTOPUS {
    tag "Octopus - $sample_id"
    // Call variants using Octopus
    // https://luntergroup.github.io/octopus/
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/octopus/", pattern: "${sample_id}.octopus.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.octopus.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.octopus.vcf.gz"), emit: vcf_file
    path "${sample_id}.octopus.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    octopus                                \
        --reference ${db}                                        \
        --reads ${bam_file}                                  \
        --regions-file ${bedfile}                                   \
        --sequence-error-model PCR.NOVASEQ             \
        --output "${sample_id}.octopus.vcf.gz"                   \
        --organism-ploidy ${params.ploidy}                        \
        --threads ${task.cpus}                                  \
        2> >(tee -a "${sample_id}.octopus.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.octopus.vcf.gz" "${sample_id}.octopus.log"
    """
}

