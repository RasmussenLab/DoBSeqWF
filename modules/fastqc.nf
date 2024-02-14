process FASTQC {
    tag "FastQC - ${sample_id}"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

    publishDir "${params.outputDir}/log/fastqc/", pattern: "${read.simpleName}_fastqc.html", mode:'copy'

    input:
    tuple val(sample_id), path(read), val(read_number)

    output:
    path "${read.simpleName}_fastqc.html", emit: fastqc_html
    path "${read.simpleName}_fastqc.zip", emit: fastqc_zip

    script:
    """
    fastqc                                  \
        --threads ${task.cpus}              \
        --quiet                             \
        ${read}
    """
    
    stub:
    """
    touch "${read.simpleName}_fastqc.html"
    touch "${read.simpleName}_fastqc.zip"
    """
}