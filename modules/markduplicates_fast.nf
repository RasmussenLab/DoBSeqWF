process MARKDUPLICATES_FAST {
    tag "Mark duplicates - $sample_id"
    // Mark duplicate reads in BAM files
    
    conda "$projectDir/envs/sambamba/environment.yaml"
    container params.container.sambamba

    publishDir "${params.outputDir}/log/markdup_sambamba/", pattern: "${sample_id}*.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val log_suffix

    output:
    tuple val(sample_id), path("${sample_id}.marked.bam"), emit: marked_bam_file
    path "${sample_id}*.log", emit: metrics_file

    script:
    def log_filename = log_suffix == "" ? "${sample_id}.markdup.log" : "${sample_id}_${log_suffix}.markdup.log"
    """
    sambamba markdup                        \
        --nthreads ${task.cpus}             \
	    --tmpdir .                          \
        --overflow-list-size 6000000        \
        ${bam_file}                         \
        ${sample_id}.marked.bam             \
        2> >(tee -a "${log_filename}" >&2)
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.markdup.log" : "${sample_id}_${log_suffix}.markdup.log"
    """
    touch "${sample_id}.marked.bam" 
    touch "${log_filename}"
    """
}