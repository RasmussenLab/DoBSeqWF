process FASTQ {
    tag "uBAM to FastQ - $sample_id"
    // Convert unaligned uBAM files to FastQ.
    
    conda "$projectDir/envs/gatk4/environment.yaml"

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_raw.fq.gz"), emit: fastq_files

    script:
    """
    gatk SamToFastq                 \
        I=${bam_file}               \
        F=${sample_id}_raw.fq.gz     \
        INTERLEAVE=true            \
        TMP_DIR=.
    """
    
    stub:
    """
    touch "${sample_id}_raw.fq.gz"
    """
}