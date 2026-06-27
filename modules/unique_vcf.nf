process UNIQUE_VCF {
    label 'process_single'
    // Merge, deduplicate, remove sample information, and index a set of VCF files.

    conda "$projectDir/envs/bcftools/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.bcftools : params.container.docker.bcftools

    publishDir "${params.outputDir}/vep_annotate/", pattern: "unique.vcf.*", mode:'copy'

    input:
    path vcf_files

    output:
    tuple path("unique.vcf.gz"), path("unique.vcf.gz.tbi"), path('variant_keys.txt'), emit: vcf_file

    script:
    """
    bcftools merge -m none --force-samples *.vcf.gz \
        | bcftools norm -d exact -O u \
        | bcftools view -G -Oz -o unique.vcf.gz
    bcftools index --tbi unique.vcf.gz
    bcftools query -f '%CHROM:%POS:%REF:%ALT\n' unique.vcf.gz > variant_keys.txt
    """

    stub:
    """
    touch "unique.vcf.gz"
    touch "unique.vcf.gz.tbi"
    touch "variant_keys.txt"
    """
}
