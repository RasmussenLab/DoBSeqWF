process MARKDUPLICATES {
    tag "Mark duplicates - $sample_id"
    // Mark duplicate reads in BAM files
    
    conda "$projectDir/envs/gatk4/environment.yaml"

    publishDir "${params.outputDir}/log/markduplicates/", pattern: "${sample_id}*.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val optical_only_tag
    val log_suffix

    output:
    tuple val(sample_id), path("${sample_id}.marked.bam"), emit: marked_bam_file
    path "${sample_id}*.log", emit: metrics_file

    script:
    def log_filename = log_suffix == "" ? "${sample_id}.dupMetric.log" : "${sample_id}_${log_suffix}.dupMetric.log"
    def optical_only = optical_only_tag ? "--duplicate-tagging-policy OpticalOnly" : ""
    """
    gatk --java-options -Xmx16g MarkDuplicates  \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        ${optical_only}                         \
	    --TMP_DIR .                             \
        -I ${bam_file}                         \
        -M ${log_filename}                     \
        -O ${sample_id}.marked.bam
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.dupMetric.log" : "${sample_id}_${log_suffix}.dupMetric.log"
    """
    touch "${sample_id}.marked.bam" 
    touch "${log_filename}"
    """
}