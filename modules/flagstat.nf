process FLAGSTAT {
    tag "Flagstat - $sample_id"
    
    conda "$projectDir/envs/samtools/environment.yaml"

    publishDir "${params.outputDir}/log/flagstat/", pattern: "${sample_id}*.flagstat", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val log_suffix

    output:
    path("${sample_id}*.flagstat"), emit: flagstat

    script:
    def log_filename = log_suffix == "" ? "${sample_id}.flagstat" : "${sample_id}_${log_suffix}.flagstat"
    """
    samtools flagstat               \
        ${bam_file}                 \
        > "${log_filename}"
    """
    
    stub:
    def log_filename = log_suffix == "" ? "${sample_id}.flagstat" : "${sample_id}_${log_suffix}.flagstat"
    """
    touch "${log_filename}"
    """
}

