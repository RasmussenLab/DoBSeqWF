process HS_METRICS {
    tag "$sample_id"
    // Collect HSMetrics for bam file

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/hs_metrics/", pattern: "${sample_id}*.hs.txt", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path(reference)
    path(bedfile)
    val log_suffix

    output:
    path("${sample_id}*.hs.txt"), emit: metrics_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def log_filename = log_suffix == "" ? "${sample_id}.hs.txt" : "${sample_id}_${log_suffix}.hs.txt"
    """
    gatk BedToIntervalList \
				 I=${bedfile} \
				 O="target_region.dict" \
				 SD=${db}

    
    gatk CollectHsMetrics                   \
        -I ${bam_file}                      \
        -O ${log_filename}                  \
        -R ${db}                            \
        --BAIT_INTERVALS "target_region.dict"        \
        --TARGET_INTERVALS "target_region.dict"      \
        --PER_TARGET_COVERAGE "target_region.dict"
    
    ## Change to correct files! Original command:
    ## --BAIT_INTERVALS {probes_intervals} \
    ## --TARGET_INTERVALS target_file_UMI_demo_data_hg38.bed \
    ## --PER_TARGET_COVERAGE probe_file_UMI_demo_data_hg38.bed
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.hs.txt" : "${sample_id}_${log_suffix}.hs.txt"
    """
    touch "${log_filename}"
    """
}