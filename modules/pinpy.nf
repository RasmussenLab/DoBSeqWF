process PINPY {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/", mode:'copy'

    input:
    path vcf_files
    path vartables
    path decodetable
    val caller

    output:
    path "pinpoint_variants"
    path "pinpoint_variants/lookup.tsv", emit: lookup_table

    script:
    """
    pin.py \
        --decodetable "${decodetable}" \
        --caller "${caller}" \
        --results-folder pinpoint_variants \
    """

    stub:
    """
    mkdir pinpoint_variants
    touch pinpoint_variants/lookup.tsv
    """
}