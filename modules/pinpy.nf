process PINPY {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/pinpoint_variants/", pattern: "results", mode:'copy'

    input:
    path vcf_files
    path vartables
    path decodetable
    path coordtable
    val caller
    val matrix_size

    output:
    path "results"
    path "results/lookup.tsv", emit: lookup_table

    script:
    """
    pin.py \
        --matrix-size "${matrix_size}" \
        --decodetable-path "${decodetable}" \
        --coordtable-path "${coordtable}" \
        --caller "${caller}" \
        --results-folder results \
    """

    stub:
    """
    mkdir results
    touch results/lookup.tsv
    """
}