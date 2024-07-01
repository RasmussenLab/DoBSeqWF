process MERGEVCFS {
    tag "MergeVCFs - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/haplotypecaller/", pattern: "${sample_id}.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}*.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}.haplotypecaller.log"

    script:
    input_files_command = vcfs.collect() {"--variant ${it}"}.join(' ')
    """
    gatk MergeVcfs  \
        ${input_files_command}      \
        -O ${sample_id}.vcf.gz \
        --TMP_DIR . \
        > >(tee -a "${sample_id}.haplotypecaller.log")  \
        2> >(tee -a "${sample_id}.haplotypecaller.log" >&2)
    """

    stub:
    input_files_command = vcfs.collect() {"--variant ${it}"}.join(' ')
    """
    echo "echo MergeVCFs - $sample_id"
    echo "gatk MergeVcfs  ${input_files_command} -O ${sample_id}.vcf.gz --TMP_DIR ."
    touch "${sample_id}.GATK.vcf.gz" "${sample_id}.haplotypecaller.log"
    """
}

