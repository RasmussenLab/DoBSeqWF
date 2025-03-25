process DEEPVARIANT {
    tag "DeepVariant - $sample_id"
    // Call variants using DeepVariant

    // Unfortunately env module on NGC is not working well. This makes the script not-portable.
    // This version works on NGC only. To run locally, comment out "singularity exec" and use "run_deepvariant" directly.

    // cpus = 8
    // memory = { 32.GB * task.attempt }
    // time = { 6.hour * task.attempt }

    publishDir "${params.outputDir}/log/deepvariant/", pattern: "${sample_id}.DV.log", mode:'copy'
    publishDir "${params.outputDir}/variants/", pattern: "${sample_id}.DV.vcf.gz", mode:'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index_file)
    path reference_genome
    path bedfile

    output:
    tuple val(sample_id), path("${sample_id}.DV.vcf.gz"), emit: vcf_file
    path "${sample_id}.DV.log"

    script:
    def db = file(params.reference_genome).getName() + ".fna"
    def ref_dir = file(params.reference_genome).getParent()
    def target_dir = file(params.bedfile).getParent()
    """
    # singularity exec --bind ${PWD} --bind ${ref_dir} --bind ${target_dir} /services/tools/deepvariant/1.5.0/deepvariant_1.5.0.sif \
    run_deepvariant \
        --ref=${db} \
        --reads=$bam_file \
        --output_vcf=${sample_id}.DV.vcf.gz \
        --model_type=WGS \
        --regions=${bedfile} \
        --intermediate_results_dir=. \
        --num_shards=${task.cpus} \
        > >(tee -a "${sample_id}.DV.log")            \
        2> >(tee -a "${sample_id}.DV.log" >&2)
    """

    stub:
    """
    touch "${sample_id}.DV.vcf.gz" "${sample_id}.DV.log"
    """
}

