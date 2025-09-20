process BGZIP {
    label 'process_low'
    tag "VCF->VCF.gz - $sample_id"
    
    conda "$projectDir/envs/samtools/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.samtools : params.container.docker.samtools

    publishDir "${params.outputDir}/norm_bgzip_idx/", pattern: "*vcf.*", mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path("${vcf_file}.gz"), path("${vcf_file}.gz.tbi"), emit: bgzf_file

    script:
    """
    bgzip ${vcf_file}
    tabix -p vcf ${vcf_file}.gz
    """
    
    stub:
    """
    touch "${vcf_file}.gz"
    touch "${vcf_file}.gz.tbi"
    """
}

