process ALIGNMENT_METRICS {
    tag "$sample_id"
    // Collect alignment metrics for bam file

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/alignment_metrics/", pattern: "${sample_id}*.align.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path(reference)
    val log_suffix

    output:
    path("${sample_id}*.align.txt"), emit: metrics_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def log_filename = log_suffix == "" ? "${sample_id}.align.txt" : "${sample_id}_${log_suffix}.align.txt"
    """
    gatk CollectAlignmentSummaryMetrics     \
        -I ${bam_file}                      \
        -O ${log_filename}                  \
        -R ${db}                            \
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.align.txt" : "${sample_id}_${log_suffix}.align.txt"
    """
    touch "${log_filename}"
    """
}