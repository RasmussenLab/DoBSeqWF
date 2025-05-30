process PILOT_PINPOINT {
    label 'process_low'
    conda "$projectDir/envs/r_env/environment.yaml"
    container params.container.r_env

    publishDir "${params.outputDir}/pinned_variants/", pattern: "*.tsv", mode:'copy'
    publishDir "${params.outputDir}/pinned_variants/outlier_plots/", pattern: "*.pdf", mode: 'copy'

    input:
    path samples, stageAs: "variants/*"
    path pooltable, stageAs: "pooltable.tsv"
    path decodetable, stageAs: "decodetable.tsv"

    output:
    path "var_id_unique_pin.tsv", emit: pinned_variants
    path "*.tsv"
    path "*.pdf"

    script:
    """
    pilot_pinpoint.R
    """

    stub:
    """
    touch var_id_unique_pin.tsv
    touch outlier_plots.pdf
    """
}