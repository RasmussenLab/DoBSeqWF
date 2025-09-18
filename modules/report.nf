process REPORT {
    label 'process_low'
    conda "$projectDir/envs/r_env/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.r_env : params.container.docker.r_env

    publishDir "${params.outputDir}/", mode:'copy', pattern: "variant_report.html"

    input:
    path rmd
    path matrix_context
    path all_variants
    path pinpoints
    path annotations
    path rescue

    output:
    path "variant_report.html"

    script:
    def annotate = annotations ? "annotations='${annotations}'," : "annotations=NULL,"
    def rescue_probs = rescue ? "rescue='${rescue}'" : "rescue=NULL"
    """
    Rscript -e \
        "rmarkdown::render(
            '${rmd}',
            params = list(
                matrix_context  = '${matrix_context}',
                pinpoints       = '${pinpoints}',
                all_variants    = '${all_variants}',
                ${annotate}
                ${rescue_probs}
            ),
        output_file = 'variant_report.html',
        quiet = TRUE
        )"
    """

    stub:
    """
    touch variant_report.html
    """
}