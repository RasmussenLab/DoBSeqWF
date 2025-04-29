process INDEX_VCF {
    label 'process_single'
    tag "Index VCF $sample_id"

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/index_vcf/", pattern: "${sample_id}.${caller}.log", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file)
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf.gz"), path("${sample_id}.${caller}.vcf.gz.tbi"), emit: vcf_w_index
    path "${sample_id}.${caller}.log"

    script:
    """
    gatk IndexFeatureFile -I ${vcf_file} \
        > >(tee -a "${sample_id}.${caller}.log")  \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf.gz.tbi"
    touch "${sample_id}.${caller}.vcf.gz"
    touch "${sample_id}.${caller}.log"
    """
}
