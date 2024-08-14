process VARTABLE_PINS {
    
    publishDir "${params.outputDir}/pinpoint_variant_tables/", mode:'copy'

    input:
    path(vcf_file)

    output:
    path("${vcf_file.baseName}.tsv"), emit: vartable

    script:
    """
    gatk VariantsToTable        \
        --split-multi-allelic   \
		--variant ${vcf_file}     \
		--output "${vcf_file.baseName}.tsv"
    """

    stub:
    """
    touch "${vcf_file.baseName}.tsv"
    """
}
