process MERGE_PINS {
    tag "Merging all pinpointables"
    // Merge pinpointables into a single VCF

    publishDir "${params.outputDir}/", mode:'copy'

    input:
    path pinpointables

    output:
    path "pinpointables.vcf"

    script:
    """
    merge_pins.sh
    """

    stub:
    """
    touch "pinpointables.vcf"
    """
}
