process INSERT_SIZE_METRICS {
    label 'process_low'
    tag "$sample_id"
    // Collect insert size metrics for bam file

    conda "$projectDir/envs/picard/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.picard : params.container.docker.picard

    publishDir "${params.outputDir}/log/insert_size_metrics/", pattern: "${sample_id}*.insert.txt", mode:'copy'
    publishDir "${params.outputDir}/log/insert_size_metrics/", pattern: "${sample_id}.pdf", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val log_suffix

    output:
    path("${sample_id}*.insert.txt"), emit: metrics_file

    script:
    def log_filename = log_suffix == "" ? "${sample_id}.insert.txt" : "${sample_id}_${log_suffix}.insert.txt"
    """
    picard CollectInsertSizeMetrics           \
        -I ${bam_file}                      \
        -O ${log_filename}                  \
        -H plot.pdf
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.insert.txt" : "${sample_id}_${log_suffix}.insert.txt"
    """
    touch "${log_filename}"
    """
}