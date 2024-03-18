process GC_METRICS {
    tag "$sample_id"
    // Collect gc bias metrics for bam file

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/gc_bias_metrics/", pattern: "${sample_id}*.gc.txt", mode:'copy'
    publishDir "${params.outputDir}/log/gc_bias_metrics/", pattern: "${sample_id}*.gc.summary.txt", mode:'copy'
    publishDir "${params.outputDir}/log/gc_bias_metrics/plot", pattern: "${sample_id}*.pdf", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path(reference)
    val log_suffix

    output:
    path("${sample_id}*.gc.txt"), emit: metrics_file
    path("${sample_id}*.gc.summary.txt"), emit: summary_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def log_filename = log_suffix == "" ? "${sample_id}.gc.txt" : "${sample_id}_${log_suffix}.gc.txt"
    def log_filename_summary = log_suffix == "" ? "${sample_id}.gc.summary.txt" : "${sample_id}_${log_suffix}.gc.summary.txt"
    """
    gatk CollectGcBiasMetrics     \
        -I ${bam_file}                      \
        -O ${log_filename}                  \
        -S ${log_filename_summary}         \
        -CHART ${sample_id}.pdf             \
        -R ${db}                            \
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.gc.txt" : "${sample_id}_${log_suffix}.gc.txt"
    def log_filename_summary = log_suffix == "" ? "${sample_id}.gc.summary.txt" : "${sample_id}_${log_suffix}.gc.summary.txt"
    """
    touch "${log_filename}"
    touch "${log_filename_summary}"
    """
}