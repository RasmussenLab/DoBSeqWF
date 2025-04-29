process MULTIQC {
    label 'process_single'
    tag "MultiQC everything"

    conda "$projectDir/envs/multiqc/environment.yaml"
    container params.container.multiqc

    publishDir "${params.outputDir}/log/multiqc/", mode:'copy'

    input:
    path(log_files)

    output:
    path "multiqc_report.html", emit: multiqc_html
    path "multiqc_data", emit: multiqc_data

    script:
    """
    multiqc . 
    """
    
    stub:
    """
    touch "multiqc_report.html"
    touch "multiqc_data"
    """
}