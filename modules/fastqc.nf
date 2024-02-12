process FASTQC {
    tag "FastQC - ${sample_id}"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

    publishDir "${params.outputDir}/log/fastqc/", pattern: "${sample_id}_fastqc.html", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc.html", emit: fastqc_html
    path "${sample_id}_fastqc.zip", emit: fastqc_zip

    script:
    """
    fastqc                                  \
        --threads ${task.cpus}              \
        --quiet                             \
        ${reads[0]} ${reads[1]}
    """
    
    stub:
    """
    touch ${sample_id}_fastqc.html
    touch ${sample_id}_fastqc.zip
    """
}