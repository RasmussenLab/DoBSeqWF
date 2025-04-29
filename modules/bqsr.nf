process BQSR {
    label 'process_low'
    tag "Calculate score recalibration - $sample_id"

    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    input:
    tuple val(sample_id), path(bam_file)
    path reference_genome
    path mills
    path g1000

    output:
    tuple val(sample_id), path(bam_file), path("${sample_id}.BQSR.table"), emit: bqsr_file

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def mills_db = file(params.mills).getName()
    def g1000_db = file(params.g1000).getName()
    """
    gatk BaseRecalibrator \
        -I ${bam_file} \
        -R ${db} \
        --known-sites ${mills_db} \
        --known-sites ${g1000_db} \
        -O ${sample_id}.BQSR.table
    """
    
    stub:
    """
    touch "${bam_file}"
    touch "${sample_id}.BQSR.table"
    """
}
