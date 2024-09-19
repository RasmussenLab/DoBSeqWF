process GENOTYPEGVCF {
    tag "GenotypeGVCFs"
    // Call variants on multisample gVCF db
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/", pattern: "genotypegvcf.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "joint.GATK.vcf.gz", mode:'copy'

    input:
    path gendb
    path reference_genome

    output:
    tuple val('joint'), path("joint.GATK.vcf.gz"), emit: vcf_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"                   \
        GenotypeGVCFs                                               \
        -R ${db}                                                    \
        -V gendb://genomicsdb                                       \
        -O "joint.GATK.vcf.gz"                                      \
        -G StandardAnnotation                                       \
        -G AS_StandardAnnotation                                    \
        -G StandardHCAnnotation                                     \
        --tmp-dir .                                                 \
        > >(tee -a "genotypegvcf.log")                              \
        2> >(tee -a "genotypegvcf.log" >&2)
    """

    stub:
    """
    touch "joint.vcf.gz" "genotypegvcf.log"
    """
}

