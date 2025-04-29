process TEST {
    label 'process_single'
    conda "$projectDir/envs/pinpy/environment.yaml"
    container params.container.pinpy

    publishDir "${params.outputDir}/test/", pattern: "test.*", mode:'copy'

    input:
    path pinned_variants
    path snv_list
    val pinpoint_method

    output:
    path "test.passed" optional true
    path "test.failed" optional true
    path "test.log"

    script:
    """
    test.py                                     \
        -t ${snv_list}                          \
        -v ${pinned_variants}                   \
        -m ${pinpoint_method}                   \
        2> >(tee -a "test.log" >&2)
    """

    stub:
    """
    touch test.passed
    touch test.log
    """
}