process BQSR {
    tag "Calculate score recalibration - $sample_id"

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

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
