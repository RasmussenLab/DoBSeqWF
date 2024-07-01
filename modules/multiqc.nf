process MULTIQC {
    tag "MultiQC everything"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

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