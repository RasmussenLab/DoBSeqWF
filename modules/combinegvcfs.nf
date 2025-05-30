process COMBINEGVCFS {
    label 'process_low'
    tag "COMBINE_GVCF"
    // Merge gVCF files into a single GenomicsDB
    
    conda "$projectDir/envs/gatk4/environment.yaml"
    container params.container.gatk

    publishDir "${params.outputDir}/log/", pattern: "genomicsdb.log", mode:'copy'
    publishDir "${params.outputDir}/splits/", pattern: "cohort.${interval}.g.vcf.gz", mode:'copy'

    input:
    path vcfs
    path index
    tuple val(interval), path(interval_file)
    path reference_genome
    path bedfile

    output:
    tuple val(interval), path("cohort.${interval}.g.vcf.gz"), path("cohort.${interval}.g.vcf.gz.tbi"), emit: gvcf_file_w_index
    path "genomicsdb.log"

    script:
    db = file(params.reference_genome).getName() + ".fna"
    input_files_command = vcfs.collect(){"--variant ${it}"}.join(' ')
    avail_mem = (task.memory.mega*0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData -XX:ConcGCThreads=1"   \
        CombineGVCFs                            \
        $input_files_command                        \
        -L ${interval_file} \
        -R ${db}      \
        -O cohort.${interval}.g.vcf.gz \
        --create-output-variant-index \
        --tmp-dir .                                                 \
        > >(tee -a "genomicsdb.log")                \
        2> >(tee -a "genomicsdb.log" >&2)
    """

    stub:
    """
    touch "genomicsdb" "genomicsdb.log"
    """
}
