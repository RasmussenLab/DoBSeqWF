process MERGEVCFS {
    tag "MergeVCFs - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    conda "$projectDir/envs/gatk4/environment.yaml"

    publishDir "${params.outputDir}/log/mergevcfs/", pattern: "${sample_id}.mergevcfs.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(vcfs)

    output:
    tuple val(sample_id), path("${sample_id}.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}_vcf.tsv", emit: vcf_info
    path "${sample_id}.mergevcfs.log"

    script:
    def input_files_command = vcfs.collect() {"--INPUT ${it}"}.join(' ')
    def publishDir = file(params.outputDir + "/variants/" + sample_id + ".GATK.vcf.gz")
    """
    gatk MergeVcfs  \
        ${input_files_command}      \
        -O ${sample_id}.GATK.vcf.gz \
        --TMP_DIR . \
        > >(tee -a "${sample_id}.mergevcfs.log")  \
        2> >(tee -a "${sample_id}.mergevcfs.log" >&2)
    
    echo -e "${sample_id}\t${publishDir}\tGATK" > "${sample_id}_vcf.tsv"
    """

    stub:
    input_files_command = vcfs.collect() {"--INPUT ${it}"}.join(' ')
    """
    echo "echo MergeVCFs - $sample_id"
    touch "${sample_id}.GATK.vcf.gz" "${sample_id}.mergevcfs.log"
    touch "${sample_id}_vcf.tsv"
    """
}

