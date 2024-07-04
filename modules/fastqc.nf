process FASTQC {
    tag "FastQC - ${sample_id}"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

    publishDir "${params.outputDir}/log/fastqc/", pattern: "${sample_id}_fastqc.html", mode:'copy'
    publishDir "${params.outputDir}/log/fastqc/", pattern: "${sample_id}_fastqc.zip", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc.html", emit: fastqc_html
    path "${sample_id}_fastqc.zip", emit: fastqc_zip

    script:
    """
    gunzip -c ${reads}                          \
        | fastqc                                \
            --threads ${task.cpus}              \
            --quiet                             \
            --noextract                         \
            stdin:"${sample_id}"
    """
    
    stub:
    """
    touch "${sample_id}_fastqc.html"
    touch "${sample_id}_fastqc.zip"
    """
}