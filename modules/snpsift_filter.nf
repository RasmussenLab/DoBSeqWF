process SNPSIFT_FILTER {
    tag "$sample_id"

    // Filter variants based on ClinVar and SNPEff annotations using SnpSift.

    publishDir "${params.outputDir}/annotated_variants/", pattern: "${sample_id}.${caller}.*.vcf", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "variants/*")
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.lof.vcf")
    tuple val(sample_id), path("${sample_id}.${caller}.pa.vcf")
    tuple val(sample_id), path("${sample_id}.${caller}.lofpa.vcf")

    script:
    """
    snpsift filter \
        "(exists LOF)" \
        ${vcf_file} \
        > ${sample_id}.${caller}.lof.vcf
    
    snpsift filter \
        "CLNSIG =~ 'Pathogenic'" \
        ${vcf_file} \
        > ${sample_id}.${caller}.pa.vcf

    snpsift filter \
        "((exists LOF) | (CLNSIG =~ 'Pathogenic'))" \
        ${vcf_file} \
        > ${sample_id}.${caller}.lofpa.vcf
    """

    stub:
    """
    touch "${sample_id}.${caller}.lof.vcf"
    touch "${sample_id}.${caller}.pa.vcf"
    touch "${sample_id}.${caller}.lofpa.vcf"
    """
}
