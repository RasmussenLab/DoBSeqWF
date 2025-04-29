process VALIDATE {
    label 'process_low'
    tag "Validate bam file - $sample_id"
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/validation/", pattern: "${sample_id}.validation.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.validation.log"), emit: validation_log

    script:
    """
    gatk ValidateSamFile                	\
        -I ${bam_file}                  	\
        -M SUMMARY                      	\
        --IGNORE_WARNINGS true          	\
        -O "${sample_id}.validation.log" 	\
				|| true
    """
    
    stub:
    """
    touch "${sample_id}.validation.log"
    """
}
