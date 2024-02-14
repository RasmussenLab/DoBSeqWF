process MARKDUPLICATES {
    tag "Mark duplicates - $sample_id"
    // Mark duplicate reads in BAM files
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/markduplicates/", pattern: "${sample_id}.dupMetric.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.marked.bam"), emit: marked_bam_file
    path "${sample_id}.dupMetric.log"

    script:
    """
    gatk MarkDuplicates                     \
        --ASSUME_SORTED true                \
        --I ${bam_file}                     \
        --O "${sample_id}.marked.bam"       \
        --M "${sample_id}.dupMetric.log"
    """
    
    stub:
    """
    touch "${sample_id}.marked.bam" 
    touch "${sample_id}.dupMetric.log"
    """
}