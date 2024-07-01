process CRISP {
    tag "CRISP - everything"
    // Call variants on all pools using CRISP
    // https://github.com/vibansal/crisp
    
    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/crisp/", pattern: "crisp.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "crisp.vcf", mode:'copy'

    input:
    path bam_files
    path reference_genome
    path bedfile

    output:
    path "crisp.vcf", emit: vcf_file
    path "crisp.log"

    script:
    def input_files_command = bam_files.collect(){"--bam ${it}"}.join(' ')
    def db = file(params.reference_genome).getName() + ".fna"
    """
    CRISP                                               \
        --ref ${db}                                     \
        ${input_files_command}                          \
        --bed ${bedfile}                                \
        --poolsize ${params.ploidy}                     \
        --VCF "crisp.vcf"                               \
        > >(tee -a "crisp.log")                         \
        2> >(tee -a "crisp.log" >&2)
    """

    stub:
    """
    touch "crisp.vcf" "crisp.log"
    """
}

