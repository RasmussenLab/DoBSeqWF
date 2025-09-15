process MATRIX_CONTEXT {
    label 'process_low'
    conda "$projectDir/envs/pinpy/environment.yaml"
    container params.container.pinpy

    publishDir "${params.outputDir}/context/", mode:'copy', pattern: "matrix_context.*"
    publishDir "${params.outputDir}/context/", mode:'copy', pattern: "decodetable.tsv"

    input:
    path pooltable
    path decodetable, stageAs: 'input_table/*'

    output:
    path "matrix_context.tsv"
    path "matrix_context.json", emit: json
    path "decodetable.tsv", emit: decodetable

    script:
    def decode = decodetable ? "--decode ${decodetable}" : ""
    """
    build_matrix.py \
        --pools ${pooltable} \
        --out-json matrix_context.json \
        --out-tsv matrix_context.tsv \
        --out-decode decodetable.tsv \
        ${decode}
    """

    stub:
    """
    touch matrix_context.json
    touch matrix_context.tsv
    touch decodetable.tsv
    """
}