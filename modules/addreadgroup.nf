process ADDREADGROUP {
    tag "Add read group - $sample_id"

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${bam_file}.bam"), emit: bam_file

    script:
    """
    gatk AddOrReplaceReadGroups     \
        --INPUT ${bam_file}         \
        --OUTPUT "${sample_id}.bam" \
        --RGID 4                    \
        --RGLB lib1                 \
        --RGPL illumina             \
        --RGPU unit1                \
        --RGSM ID ${sample_id}
    """
    
    stub:
    """
    touch "${bam_file}.bam"
    """
}

