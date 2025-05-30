process HAPLOTYPECALLER_SPLIT {
    label 'process_low'
    tag "HaplotypeCaller - $sample_id"
    // Call variants using GATK - HaplotypeCaller
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/haplotypecaller/${sample_id}/", pattern: "${sample_id}*.g.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/variants/gvcf/${sample_id}/", pattern: "${sample_id}*.GATK.g.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file), val(interval)
    path reference_genome
    path interval_list

    output:
    tuple val(sample_id), val(interval), path("${sample_id}*.GATK.g.vcf.gz"), path("${sample_id}*.GATK.g.vcf.gz.tbi"), emit: gvcf_file
    path "${sample_id}*.g.haplotypecaller.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def id = interval ? "${sample_id}.${interval}" : sample_id
    def add_intervals = interval ? "-L ${interval}" : ""
    def add_intersection = interval ? "--interval-set-rule INTERSECTION" : ""
    """
    gatk --java-options "-Xmx8g -XX:-UsePerfData" HaplotypeCaller   \
        -R ${db}                                                    \
        -I ${bam_file}                                              \
        -L ${interval_list}                             \
        ${add_intervals}                                \
        ${add_intersection}                             \
        -O "${id}.GATK.g.vcf.gz"                             \
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
        --max-reads-per-alignment-start 0                           \
        --create-output-variant-index \
        --tmp-dir .                                                 \
        > >(tee -a "${id}.g.haplotypecaller.log")            \
        2> >(tee -a "${id}.g.haplotypecaller.log" >&2)
    """

    stub:
    def id = interval ? "${sample_id}.${interval}" : sample_id
    """
    touch "${id}.GATK.g.vcf.gz" "${id}.GATK.g.vcf.gz.tbi" "${id}.g.haplotypecaller.log"
    """
}

