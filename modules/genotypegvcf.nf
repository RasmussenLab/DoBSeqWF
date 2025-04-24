process GENOTYPEGVCF {
    tag "GenotypeGVCFs"
    // Call variants on multisample gVCF db
    
    conda "$projectDir/envs/gatk4/environment.yaml"

    publishDir "${params.outputDir}/log/genotypegvcf/", pattern: "${sample_id}.genotypegvcf.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(gvcf_file), path(index_file)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}.GATK.vcf.gz"), emit: vcf_file
    path "${sample_id}_vcf.tsv", emit: vcf_info
    path "${sample_id}.genotypegvcf.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def publishDir = file(params.outputDir + "/variants/" + sample_id + ".GATK.vcf.gz")
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"                   \
        GenotypeGVCFs                                               \
        -R ${db}                                                    \
        -V ${gvcf_file}                                             \
        -O "${sample_id}.GATK.vcf.gz"                               \
        -G StandardAnnotation                                       \
        -G AS_StandardAnnotation                                    \
        -G StandardHCAnnotation                                     \
        --sample-ploidy ${params.ploidy}                            \
        --tmp-dir .                                                 \
        > >(tee -a "${sample_id}.genotypegvcf.log")                              \
        2> >(tee -a "${sample_id}.genotypegvcf.log" >&2)    
    
    echo -e "${sample_id}\t${publishDir}\tGATK" > "${sample_id}_vcf.tsv"
    """

    stub:
    """
    touch "${sample_id}.GATK.vcf.gz" "${sample_id}.genotypegvcf.log"
    touch "${sample_id}_vcf.tsv"
    """
}

