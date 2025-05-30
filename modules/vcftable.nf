process VCFTABLE {
    label 'process_single'
    publishDir "${projectDir}/", pattern: "vcftable.tsv", mode:'copy'

    input:
    val(vcf_info)

    output:
    path("vcftable.tsv")

    script:
    """
    cat ${vcf_info.join(' ').trim()} > vcftable.tsv
    """

    stub:
    """
    touch vcftable.tsv
    """
}