process SNPEFF {
    tag "$sample_id"

    // Add functional annotations to a VCF file using SnpEff.
    // -strict and -canon options ensure that only high-confidence annotations are added.
    // -lof option to add loss-of-function annotations.
    // -noMotif avoids adding motif annotations.

    conda "$projectDir/envs/snpeff/environment.yaml"

    publishDir "${params.outputDir}/log/snpeff/${sample_id}/", pattern: "${sample_id}.${caller}.log", mode:'copy'
    publishDir "${params.outputDir}/log/snpeff/stats/", pattern: "${sample_id}.${caller}_snpeff_stats.csv", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "variants/*")
    val caller
    val snpeff_db
    path snpeff_cache
    path snpeff_config

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: snpeff_vcf
    path "${sample_id}.${caller}_snpeff_stats.csv", emit: snpeff_stats
    path "${sample_id}.${caller}.log"

    script:
    """
    ${params.snpeff}                                        \
        ${snpeff_db}                                        \
        -c ${snpeff_config}                                 \
        -csvStats "${sample_id}.${caller}_snpeff_stats.csv" \
        -lof                                                \
        -strict                                             \
        -noMotif                                            \
        -canon                                              \
        -dataDir ${snpeff_cache}                            \
        ${vcf_file}                                         \
        > ${sample_id}.${caller}.vcf              \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    touch "${sample_id}.${caller}_snpeff_stats.csv"
    touch "${sample_id}.${caller}.log"
    """
}
