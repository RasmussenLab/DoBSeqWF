process ALIGNMENT {
    tag "Alignment - $sample_id"
    // Align reads to the reference genome using BWA, convert to BAM and sort.

    conda "$projectDir/envs/bwa/environment.yaml"

    publishDir "${params.outputDir}/log/mapping/", pattern: "${sample_id}.log", mode:'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: raw_bam_file
    path "${sample_id}.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    // int bwa_cpus = Math.max((task.cpus * 0.75) as int, 1)
    // int samtools_cpus = Math.max(task.cpus - bwa_cpus, 1)

    """
    bwa-mem2 mem                                                        \
        -t ${task.cpus-10}                                              \
        -K 100000000                                                    \
        -v 2                                                            \
        ${db}                                                           \
        ${reads[0]}                                                     \
        ${reads[1]}                                                     \
        2> >(tee -a "${sample_id}.log" >&2)                             \
            | samtools sort -@ 10 -o "${sample_id}.bam" -
    """
    
    stub:
    """
    touch "${sample_id}.bam"
    touch "${sample_id}.log"
    """
}
