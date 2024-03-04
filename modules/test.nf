process TEST {

    // cpus = { 2 * task.attempt }
    // memory = { 4.GB * task.attempt }
    // time = { 1.hour * task.attempt }

    publishDir "${params.outputDir}/test/", pattern: "test.*", mode:'copy'

    input:
    path pinned_variants
    path snv_list

    output:
    path "test.passed" optional true
    path "test.failed" optional true
    path "test.log"

    script:
    """
    test.py                                     \
        -t ${snv_list}                          \
        -v ${pinned_variants}                   \
        2> >(tee -a "test.log" >&2)
    """

    stub:
    """
    touch test.passed
    touch test.log
    """
}