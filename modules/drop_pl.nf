process DROP_PL {
    label 'process_single'
    tag "$sample_id"
    // Drop PL tag from VCF file.

    conda "$projectDir/envs/bcftools/environment.yaml"
    container params.container.bcftools

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "input/*"), path(index, stageAs: "input/*")
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: vcf_file

    script:
    """
    bcftools annotate                   \
        -x "FORMAT/PL"                  \
        -Ov                             \
        -o ${sample_id}.${caller}.vcf   \
        ${vcf_file}
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    """
}
