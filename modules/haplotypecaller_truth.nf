process HC_TRUTH {
    tag "HaplotypeCaller - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/haplotypecaller/", pattern: "${sample_id}.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"       \
        HaplotypeCaller                                 \
        -R ${db}                                        \
        -I ${bam_file}                                  \
        -L ${bedfile}                                   \
        -O "${sample_id}.GATK.vcf.gz"                   \
        > >(tee -a "${sample_id}.haplotypecaller.log")  \
        2> >(tee -a "${sample_id}.haplotypecaller.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.GATK.vcf.gz" "${sample_id}.haplotypecaller.log"
    """
}

