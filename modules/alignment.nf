process ALIGNMENT {
    tag "Alignment - $sample_id"
    // Align reads to the reference genome using BWA, convert to BAM and sort.
    // cpus = { 15 * task.attempt }
    // memory = { 56.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/mapping/", pattern: "${sample_id}.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}.raw.bam"), emit: raw_bam_file
    path "${sample_id}.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    bwa-mem2 mem                                                        \
        -t ${task.cpus}                                                 \
        -K 100000000                                                    \
        -v 2                                                            \
        ${db}                                                           \
        ${reads[0]}                                                     \
        ${reads[1]}                                                     \
        2> >(tee -a "${sample_id}.log" >&2)                             \
            | samtools sort -o "${sample_id}.raw.bam" -
    """
    
    stub:
    """
    touch "${sample_id}.raw.bam"
    touch "${sample_id}.log"
    """
}