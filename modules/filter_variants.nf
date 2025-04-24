process FILTER_VARIANTS {

    conda "$projectDir/envs/filter/environment.yaml"

    publishDir "${params.outputDir}/filtered_variants/", pattern: "*.vcf", mode:'copy'

    input:
    path vcf_files, stageAs: "input/*"
    path snv_model
    path indel_model
    path model_class, stageAs: "models/*"
    val filter_method
    val filter_indels

    output:
    path("*.vcf"), emit: filtered_vcfs

    script:
    def input_vcfs = vcf_files.join(' ')
    def filter = filter_indels ? "" : "--indel-threshold 1"
    """
    export PYTHONPATH=\${PYTHONPATH:+\$PYTHONPATH:}\$PWD
    dwf_filter.py                           \
        --input ${input_vcfs}               \
        --output ./                         \
        --snv-model ${snv_model}            \
        --indel-model ${indel_model}        \
        --threshold-type ${filter_method}   \
        ${filter}
    """

    stub:
    """
    for vcf in ${vcf_files}; do
        touch \$(basename \$vcf)
    done
    """
}
