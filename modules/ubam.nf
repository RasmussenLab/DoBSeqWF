process UBAM {
    tag "FastQ to uBAM - $sample_id"
    // Convert FastQ files to unaligned BAM files.
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.ubam"), emit: unaligned_bam_file

    script:
    """
    gatk FastqToSam             \
        O="${sample_id}.ubam"   \
        F1=${reads[0]}          \
        F2=${reads[1]}          \
        SM=${sample_id}         \
        LB=Library1             \
        PU=Unit1                \
        PL=Illumina             \
        TMP_DIR=.
    """
    
    stub:
    """
    touch "${sample_id}.ubam"
    """
}