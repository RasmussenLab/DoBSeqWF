process GENOTYPEGVCF_SPLIT {
    tag "GenotypeGVCFs"
    // Call variants on multisample gVCF db
    
    conda "$projectDir/envs/gatk4/environment.yaml"

    publishDir "${params.outputDir}/log/genotypegvcf/${sample_id}/", pattern: "${sample_id}*.genotypegvcf.log", mode:'copy'

    input:
    tuple val(sample_id), val(interval), path(vcf_file), path(vcf_index)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}*.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}*.genotypegvcf.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def id = interval ? "${sample_id}.${interval}" : sample_id
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"                   \
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

