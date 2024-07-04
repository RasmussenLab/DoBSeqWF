process NORMALISE_VCF {
    tag "$sample_id"
    // 1. Splits multi-allelic variants into multiple lines
    // 2. Left-aligns indels (standard for GATK).
    // 3. Checks refs against reference genome

    publishDir "${params.outputDir}/log/normalise_vcf/${sample_id}/", pattern: "${sample_id}.${caller}.log", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file), path(index)
    path reference_genome
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: norm_vcf
    path "${sample_id}.${caller}.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    bcftools norm                           \
        --check-ref e                       \
        --fasta-ref "${db}"                 \
        -m -both                            \
        -Ov                                 \
        -o "${sample_id}.${caller}.vcf"     \
        "${vcf_file}"                       \
        > >(tee -a "${sample_id}.${caller}.log")  \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    touch "${sample_id}.${caller}.log"
    """
}
