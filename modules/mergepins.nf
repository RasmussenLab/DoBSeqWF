process MERGE_PINS {
    tag "Merging all pinpointables"
    // Merge pinpointables into a single VCF

    conda "$projectDir/envs/bcftools/environment.yaml"
    container params.container.bcftools

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
