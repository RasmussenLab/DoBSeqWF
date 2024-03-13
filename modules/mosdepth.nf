process MOSDEPTH {
    tag "mosdepth - ${sample_id}"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.per-base.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.regions.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.mosdepth.global.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.mosdepth.region.dist.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path targets_bed
    val log_suffix

    output:
    path "${sample_id}*.per-base.bed.gz", emit: per_base_depth
    path "${sample_id}*.regions.bed.gz", emit: region_depth
    path "${sample_id}*.mosdepth.global.dist.txt", emit: global_dist
    path "${sample_id}*.mosdepth.region.dist.txt", emit: region_dist

    script:
    def log_filename = log_suffix == "" ? "${sample_id}" : "${sample_id}_${log_suffix}"
    """
    mosdepth                                \
        --threads ${task.cpus}              \
        --by ${targets_bed}                 \
        ${log_filename}                     \
        ${bam}
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}" : "${sample_id}_${log_suffix}"
    """
    touch ${log_filename}.per-base.bed.gz
    touch ${log_filename}.regions.bed.gz
    touch ${log_filename}.mosdepth.global.dist.txt
    touch ${log_filename}.mosdepth.region.dist.txt
    """
}
