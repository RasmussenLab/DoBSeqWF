process HAPLOTYPECALLER {
    label 'process_low'
    tag "HaplotypeCaller - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/haplotypecaller/", pattern: "${sample_id}.g.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.GATK.g.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.GATK.g.vcf.gz"), path("${sample_id}.GATK.g.vcf.gz.tbi"), emit: gvcf_file
    path "${sample_id}.g.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData" HaplotypeCaller   \
        -R ${db}                                                    \
        -I ${bam_file}                                              \
        -L ${bedfile}                                               \
        -O "${sample_id}.GATK.g.vcf.gz"                             \
        -ploidy ${params.ploidy}                                    \
        -ERC GVCF                                                   \
        -G StandardAnnotation                                       \
        -G AS_StandardAnnotation                                    \
        -G StandardHCAnnotation                                     \
        -A GcContent                                                \
        -A HmerIndelLength                                          \
        -A IndelLength                                              \
        -A AssemblyComplexity                                       \
        -A LikelihoodRankSumTest                                    \
        -A ClippingRankSumTest                                      \
        -A HaplotypeFilteringAnnotation                             \
        -A GenotypeSummaries                                        \
        -A AlleleFraction                                           \
        --disable-read-filter NotDuplicateReadFilter                \
        --max-alternate-alleles 3                                   \
        --max-num-haplotypes-in-population 1000                     \
        --max-reads-per-alignment-start 0                           \
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

