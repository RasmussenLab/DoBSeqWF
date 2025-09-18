process PIN_BASIC {
    label 'process_low'
    conda "$projectDir/envs/filter_variants/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.filter_variants : params.container.docker.filter_variants

    publishDir "${params.outputDir}/variant_compilation/", mode:'copy', pattern: "*variants.tsv"

    input:
    path vcf_files
    path matrix_context

    output:
    path "all_pool_variants.tsv", emit: pool_variants
    path "pinpoint_variants.tsv", emit: pinpoint_variants

    script:
    """
    pin_basic.py                    \
        --context ${matrix_context} \
        --vcf-folder .              \
        --output .
    """

    stub:
    """
    touch all_pool_variants.tsv
    touch pinpoint_variants.tsv
    """
}