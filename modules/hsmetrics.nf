process HSMETRICS {
    tag "HSMETRICS - $sample_id"
    // Collect HSMetrics for bam file

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    input:
    tuple val(sample_id), path(bam_file)
    path(reference)
    path(bed_file)

    output:
    path("${sample_id}.txt"), emit: metrics_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk CollectHsMetrics                   \
        -I ${bam_file}                      \
        -O ${sample_id}.txt                 \
        -R ${db}                            \
        --BAIT_INTERVALS ${bed_file}        \
        --TARGET_INTERVALS ${bed_file}      \
        --PER_TARGET_COVERAGE ${bed_file}
    
    ## Change to correct files! Original command:
    ## --BAIT_INTERVALS {probes_intervals} \
    ## --TARGET_INTERVALS target_file_UMI_demo_data_hg38.bed \
    ## --PER_TARGET_COVERAGE probe_file_UMI_demo_data_hg38.bed
    """
    
    stub:
    """
    touch "${sample_id}.ubam"
    """
}