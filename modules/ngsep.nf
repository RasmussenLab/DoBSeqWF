process NGSEP {
    label 'process_low'
    tag "NGSEP - $sample_id"
    // Call variants using NGSEP

    conda "$projectDir/envs/ngsep/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.ngsep : params.container.docker.ngsep

    publishDir "${params.outputDir}/log/ngsep/", pattern: "${sample_id}.ngsep.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.ngsep.vcf", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.ngsep.vcf"), emit: vcf_file
    path "${sample_id}.ngsep.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    ngsep \
        SingleSampleVariantsDetector \
        -maxAlnsPerStartPos 50000 \
        -h 0.1 \
        -ploidy ${params.ploidy} \
        -r ${db} \
        -i ${bam_file} \
        -o "${sample_id}.ngsep" \
        > >(tee -a "${sample_id}.ngsep.log")            \
        2> >(tee -a "${sample_id}.ngsep.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.ngsep.vcf" "${sample_id}.ngsep.log"
    """
}

