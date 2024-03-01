process CRAMTABLE {

    publishDir "${projectDir}/", pattern: "cramtable.tsv", mode:'copy'

    input:
    val(cram_info)

    output:
    path("cramtable.tsv")

    script:
    """
    cat ${cram_info.join(' ').trim()} > cramtable.tsv
    """

    stub:
    """
    touch cramtable.tsv
    """
}