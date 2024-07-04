process VARTABLE {
    tag "VCF to TSV $sample_id"

    publishDir "${params.outputDir}/log/variant_tables/${sample_id}/", pattern: "${sample_id}.${caller}.log", mode:'copy'
    publishDir "${params.outputDir}/variant_tables/", pattern: "${sample_id}.${caller}.tsv", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file), path(index)
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.tsv"), emit: vartable
    path "${sample_id}.${caller}.log"

    script:
    """
    gatk VariantsToTable        \
        --split-multi-allelic   \
		--variant ${vcf_file}     \
		--output "${sample_id}.${caller}.tsv" \
        > >(tee -a "${sample_id}.${caller}.log")  \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.tsv"
    touch "${sample_id}.${caller}.log"
    """
}
