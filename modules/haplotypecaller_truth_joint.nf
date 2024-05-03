process HC_TRUTH_JOINT {
    tag "HaplotypeCaller - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/haplotypecaller/", pattern: "${sample_id}.g.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.g.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    path "${sample_id}.GATK.g.vcf.gz", emit: gvcf_file
    path "${sample_id}.GATK.g.vcf.gz.tbi", emit: gvcf_index
    path "${sample_id}.g.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData"                   \
        HaplotypeCaller                                             \
        -R ${db}                                                    \
        -L ${bedfile}                                               \
        -I ${bam_file}                                              \
        -O "${sample_id}.GATK.g.vcf.gz"                             \
        -ERC GVCF                                                   \
        -G StandardAnnotation                                       \
        -G AS_StandardAnnotation                                    \
        -G StandardHCAnnotation                                     \
        --create-output-variant-index                               \
        --tmp-dir .                                                 \
        > >(tee -a "${sample_id}.g.haplotypecaller.log")            \
        2> >(tee -a "${sample_id}.g.haplotypecaller.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.GATK.g.vcf.gz" "${sample_id}.g.haplotypecaller.log"
    """
}

