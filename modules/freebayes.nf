process FREEBAYES {
    label 'process_low'
    tag "Freebayes call - $sample_id"
    // Call variants using Lofreq - call parallel
    
    conda "$projectDir/envs/freebayes/environment.yaml"
    container params.container.freebayes

    publishDir "${params.outputDir}/log/freebayes/", pattern: "${sample_id}.freebayes.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.freebayes.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.freebayes.vcf.gz"), emit: vcf_file
    path "${sample_id}_vcf.tsv", emit: vcf_info
    path "${sample_id}.freebayes.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def publishDir = file(params.outputDir + "/variants/" + sample_id + ".freebayes.vcf.gz")
    """
    freebayes                                       \
        --ploidy ${params.ploidy}                   \
        --use-best-n-alleles 3                      \
        --min-alternate-fraction 0.005              \
        --pooled-discrete                           \
        --min-alternate-count 2                     \
        --report-all-haplotype-alleles              \
        -f ${db}                                    \
        -t ${bedfile}                               \
        -v "${sample_id}.freebayes.vcf.gz"          \
        ${bam_file}                                 \
        > >(tee -a "${sample_id}.freebayes.log")       \
        2> >(tee -a "${sample_id}.freebayes.log" >&2)

    echo -e "${sample_id}\t${publishDir}\tfreebayes" > "${sample_id}_vcf.tsv"
    """
    
    stub:
    """
    touch "${sample_id}.freebayes.vcf.gz"
    touch "${sample_id}.freebayes.log"
    touch "${sample_id}_vcf.tsv"
    """
}

