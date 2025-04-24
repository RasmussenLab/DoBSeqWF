process SNPSIFT_CLINVAR {
    tag "$sample_id"

    // Add ClinVar annotations to a VCF file using SnpSift.

    conda "$projectDir/envs/snpeff/environment.yaml"

    publishDir "${params.outputDir}/log/snpsift_clinvar/${sample_id}/", pattern: "${sample_id}.${caller}.log", mode:'copy'
    publishDir "${params.outputDir}/annotated_variants/", pattern: "${sample_id}.${caller}.vcf", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "variants/*")
    val caller
    path clinvar_db

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: clinvar_vcf
    path "${sample_id}.${caller}.log"

    script:
    def db = file(params.clinvar_db).getName()
    """
    ${params.snpsift} annotate                          \
        ${db}                                           \
        ${vcf_file}                                     \
        > ${sample_id}.${caller}.vcf          \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    touch "${sample_id}.${caller}.log"
    """
}
