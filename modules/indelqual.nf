process INDELQUAL {
    tag "Lofreq indelqual - $sample_id"
    // Tag indel quality using Lofreq - indelqual
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/lofreq_indelqual/", pattern: "${sample_id}.lofreq.iq.log", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path reference_genome

    output:
    tuple val(sample_id), path("${sample_id}.iq.bam"), emit: iq_bam_file
    path "${sample_id}.lofreq.iq.log"

    script:
    """
    lofreq indelqual                        \
        --dindel                            \
        -f ${reference_genome}              \
        -o "${sample_id}.iq.bam"            \
        ${bam_file}                         \
        &> "${sample_id}.lofreq.iq.log"
    """
    
    stub:
    """
    touch "${sample_id}.iq.bam" "${sample_id}.lofreq.iq.log"
    """
}

