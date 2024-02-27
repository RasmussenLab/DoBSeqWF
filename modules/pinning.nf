process PINNING {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/pinned_variants/", pattern: "*.tsv", mode:'copy'
    publishDir "${params.outputDir}/pinned_variants/outlier_plots/", pattern: "*.pdf", mode: 'copy'

    input:
    path samples, stageAs: "variants/*"
    path pooltable, stageAs: "pooltable.tsv"
    path decodetable, stageAs: "decodetable.tsv"

    output:
    path "*.tsv", emit: pinned_variants
    path "*.pdf"

    script:
    """
    pinning.R \
        --nextflow
    """

    stub:
    """
    touch var_id_unique_pin.tsv
    touch outlier_plots.pdf
    """
}