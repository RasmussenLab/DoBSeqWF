process SNPSIFT_FILTER {
    label 'process_low'
    tag "$sample_id"

    // Filter variants based on ClinVar and SNPEff annotations using SnpSift.

    conda "$projectDir/envs/snpeff/environment.yaml"
    container params.container.snpeff

    publishDir "${params.outputDir}/annotated_variants/lof/", pattern: "${sample_id}.${caller}.lof.vcf", mode:'copy'
    publishDir "${params.outputDir}/annotated_variants/pathogenic/", pattern: "${sample_id}.${caller}.p.vcf", mode:'copy'
    publishDir "${params.outputDir}/annotated_variants/lofp/", pattern: "${sample_id}.${caller}.lofp.vcf", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "variants/*")
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.lof.vcf")
    tuple val(sample_id), path("${sample_id}.${caller}.p.vcf")
    tuple val(sample_id), path("${sample_id}.${caller}.lofp.vcf")

    script:
    """
    ${params.snpsift} filter \
        "(exists LOF)" \
        ${vcf_file} \
        > ${sample_id}.${caller}.lof.vcf
    
    ${params.snpsift} filter \
        "CLNSIG =~ 'Pathogenic'" \
        ${vcf_file} \
        > ${sample_id}.${caller}.p.vcf

    ${params.snpsift} filter \
        "((exists LOF) | (CLNSIG =~ 'Pathogenic'))" \
        ${vcf_file} \
        > ${sample_id}.${caller}.lofp.vcf
    """

    stub:
    """
    touch "${sample_id}.${caller}.lof.vcf"
    touch "${sample_id}.${caller}.p.vcf"
    touch "${sample_id}.${caller}.lofp.vcf"
    """
}
