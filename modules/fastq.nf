process FASTQ {
    tag "uBAM to FastQ - $sample_id"
    // Convert unaligned uBAM files to FastQ.
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}.fq.gz"), emit: fastq_files

    script:
    """
    gatk SamToFastq                 \
        I=${bam_file}               \
        F=${sample_id}_R1.fq.gz     \
        F2=${sample_id}_R2.fq.gz    \
    """
    
    stub:
    """
    touch "${sample_id}_R1.fq.gz"
    touch "${sample_id}_R2.fq.gz"
    """
}