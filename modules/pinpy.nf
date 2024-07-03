process PINPY {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/pinpoint_variants/", mode:'copy'

    input:
    path vcf_files
    path vartables
    path decodetable
    val caller

    output:
    path "results"
    path "results/lookup.tsv", emit: lookup_table

    script:
    """
    pin.py \
        --decodetable "${decodetable}" \
        --caller "${caller}" \
        --results-folder results \
    """

    stub:
    """
    mkdir results
    touch results/lookup.tsv
    """
}