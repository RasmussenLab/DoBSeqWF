process DISCARD {
    label 'process_single'
    publishDir "${params.outputDir}/undiscarded_variants/", pattern: "*.vcf", mode:'copy'

    input:
    path vcf_file, stageAs: "input/*"

    output:
    path("*.vcf"), emit: vcf_file

    shell:
    """
    awk '/^#/||\$7!="ML_FILTERED"' $vcf_file > \$(basename $vcf_file)
    """

    stub:
    """
    touch \$(basename $vcf_file)
    """
}
