process FASTQC {
    tag "FastQC - ${sample_id}"

    conda "$projectDir/envs/fastqc/environment.yaml"
    container params.container.fastqc

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