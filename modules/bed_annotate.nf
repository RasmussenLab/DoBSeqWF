process BED_ANNOTATE {
    tag "$sample_id"
    // Annotate VCF by third column in BED file

    publishDir "${params.outputDir}/log/bed_annotate/${sample_id}/", pattern: "${sample_id}.${caller}.log", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "input/*")
    path bedfile
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: vcf_file
    path "${sample_id}.${caller}.log"

    script:
    """
    bcftools annotate \
        -a ${bedfile} \
        -c CHROM,FROM,TO,GENE \
        -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name from target bedfile">') \
        -o ${sample_id}.${caller}.vcf \
        ${vcf_file} \
        > >(tee -a "${sample_id}.${caller}.log")  \
        2> >(tee -a "${sample_id}.${caller}.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    touch "${sample_id}.${caller}.log"
    """
}
