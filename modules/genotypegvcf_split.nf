process GENOTYPEGVCF_SPLIT {
    label 'process_low'
    tag "GenotypeGVCFs"
    // Call variants on multisample gVCF db
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.gatk : params.container.docker.gatk

    publishDir "${params.outputDir}/log/genotypegvcf/${sample_id}/", pattern: "${sample_id}*.genotypegvcf.log", mode:'copy'

    input:
    tuple val(sample_id), val(interval), path(vcf_file), path(vcf_index)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}*.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}*.genotypegvcf.log"

    script:
    def db = file(params.reference_genome).name
    def id = interval ? "${sample_id}.${interval}" : sample_id
    def avail_mem = (task.memory.mega*0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData"                   \
        GenotypeGVCFs                                               \
        -R ${db}                                                    \
        -V ${vcf_file}                                              \
        -O "${id}.GATK.vcf.gz"                                      \
        -stand-call-conf 0                                          \
        -G StandardAnnotation                                       \
        -G AS_StandardAnnotation                                    \
        -G StandardHCAnnotation                                     \
        --sample-ploidy ${params.ploidy}                            \
        --tmp-dir .                                                 \
        > >(tee -a "${id}.genotypegvcf.log")                        \
        2> >(tee -a "${id}.genotypegvcf.log" >&2)
    """

    stub:
    def id = interval ? "${sample_id}.${interval}" : sample_id
    """
    touch "${id}.GATK.vcf.gz" "${id}.genotypegvcf.log"
    """
}

