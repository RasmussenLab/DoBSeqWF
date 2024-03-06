process MOSDEPTH {
    tag "mosdepth - ${sample_id}"

    // cpus = { 16 * task.attempt }
    // memory = { 10.GB * task.attempt }
    // time = { 4.hour * task.attempt }

    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}.per-base.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}.regions.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}.mosdepth.global.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/log/mosdepth/", pattern: "${sample_id}.mosdepth.region.dist.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path targets_bed

    output:
    path "${sample_id}.per-base.bed.gz", emit: per_base_depth
    path "${sample_id}.regions.bed.gz", emit: region_depth
    path "${sample_id}.mosdepth.global.dist.txt", emit: global_dist
    path "${sample_id}.mosdepth.region.dist.txt", emit: region_dist

    script:
    """
    mosdepth                                \
        --threads ${task.cpus}              \
        --by ${targets_bed}                 \
        ${sample_id}                        \
        ${bam}
    """
    
    stub:
    """
    touch ${sample_id}.per-base.bed.gz
    touch ${sample_id}.regions.bed.gz
    touch ${sample_id}.mosdepth.global.dist.txt
    touch ${sample_id}.mosdepth.region.dist.txt
    """
}
