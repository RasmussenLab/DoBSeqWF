process DUPLICATE_METRICS {
    tag "$sample_id"
    // Mark duplicate reads in BAM files
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/duplicate_metrics/", pattern: "${sample_id}*.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val log_suffix

    output:
    path "${sample_id}*.txt", emit: metrics_file

    script:
    def log_filename = log_suffix == "" ? "${sample_id}.dup.txt" : "${sample_id}_${log_suffix}.dup.txt"
    """
    gatk MarkDuplicates                     \
        --I ${bam_file}                     \
        --M ${log_filename}                 \
        --O ${sample_id}.marked.bam
    
    rm ${sample_id}.marked.bam
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.dup.txt" : "${sample_id}_${log_suffix}.dup.txt"
    """
    touch "${log_filename}"
    """
}