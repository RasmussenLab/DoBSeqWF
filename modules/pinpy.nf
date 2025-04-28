process PINPY {

    conda "$projectDir/envs/pinpy/environment.yaml"
    container params.container.pinpy

    publishDir "${params.outputDir}/", mode:'copy'

    input:
    path vcf_files
    path vartables
    path decodetable
    val caller

    output:
    path "pinpoint_variants"
    path "pinpoint_variants/lookup.tsv", emit: lookup_table
    path "pinpoint_variants/all_pins/*.vcf", emit: vcf_all_pins
    path "pinpoint_variants/unique_pins/*.vcf", emit: vcf_unique_pins
    path "pinpoint_variants/unique_pins/*.vcf.gz", emit: vcf_unique_2d_pins

    script:
    """
    pin.py \
        --decodetable "${decodetable}" \
        --caller "${caller}" \
        --results-folder pinpoint_variants \
    """

    stub:
    """
    mkdir -p pinpoint_variants/all_pins/ pinpoint_variants/unique_pins/
    touch pinpoint_variants/lookup.tsv
    for vcf in ${vcf_files}; do
        touch pinpoint_variants/all_pins/\$(basename \$vcf)
        touch pinpoint_variants/unique_pins/\$(basename \$vcf)
        touch pinpoint_variants/unique_pins/\$(basename \$vcf).gz
    done
    """
}