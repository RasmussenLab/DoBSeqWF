process MOSDEPTH {
    label 'process_low'
    tag "mosdepth - ${sample_id}"

    // Added additional high coverage thresholds:
    // 5x allele depth (240), 10x allele depth (480), 30x allele depth (1440), 50x allele depth (2400)

    conda "$projectDir/envs/mosdepth/environment.yaml"
    container params.container.mosdepth

    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.per-base.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.regions.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.mosdepth.global.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.mosdepth.region.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}*.mosdepth.summary.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path targets_bed
    val log_suffix

    output:
    path "${sample_id}*.per-base.bed.gz", emit: per_base_depth
    path "${sample_id}*.regions.bed.gz", emit: region_depth
    path "${sample_id}*.mosdepth.global.dist.txt", emit: global_dist
    path "${sample_id}*.mosdepth.region.dist.txt", emit: region_dist
    path "${sample_id}*.mosdepth.summary.txt", emit: summary

    script:
    def log_filename = log_suffix == "" ? "${sample_id}" : "${sample_id}_${log_suffix}"
    """
    mosdepth                                \
        --threads ${task.cpus}              \
        --by ${targets_bed}                 \
        --thresholds 1,10,30,240,480,1440,2400    \
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
    touch ${log_filename}.mosdepth.summary.txt
    """
}
