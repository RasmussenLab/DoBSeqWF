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
    tuple val(sample_id), path("${sample_id}.iq.bam"), emit: bam_file
    path "${sample_id}.lofreq.iq.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    """
    lofreq indelqual                                    \
        --dindel                                        \
        -f ${db}                                        \
        -o "${sample_id}.iq.bam"                        \
        ${bam_file}                                     \
        > >(tee -a "${sample_id}.lofreq.iq.log")        \
        2> >(tee -a "${sample_id}.lofreq.iq.log" >&2)
    """
    
    stub:
    """
    touch "${sample_id}.iq.bam"
    touch "${sample_id}.lofreq.iq.log"
    """
}

