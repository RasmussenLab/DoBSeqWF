process RESCUE {
    label 'process_low'
    conda "$projectDir/envs/rescue/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.marbl : params.container.docker.marbl


    publishDir "${params.outputDir}/rescue_probabilities/", mode:'copy', pattern: "predictions.tsv"

    input:
    path sampletable
    path decodetable
    path vcf_files
    path mpileup

    output:
    path "predictions.tsv", emit: probabilities

    script:
    def decode = decodetable ? "--decodetable ${decodetable}" : ""
    """
    marbl                               \
        --vcf-folder .                  \
        --mpileup-folder .              \
        --sampletable ${sampletable}    \
        ${decode}                       \
        --output .
    """

    stub:
    """
    touch predictions.tsv
    """
}
