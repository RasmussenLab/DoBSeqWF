process DROP_PL {
    label 'process_single'
    tag "$sample_id"
    // Drop PL tag from VCF file.

    conda "$projectDir/envs/bcftools/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.bcftools : params.container.docker.bcftools

    input:
    tuple val(sample_id), path(vcf_file, stageAs: "input/*"), path(index, stageAs: "input/*")
    val caller

    output:
    tuple val(sample_id), path("${sample_id}.${caller}.vcf"), emit: vcf_file

    script:
    """
    # 1) Drop PL from the VCF
    bcftools annotate \\
        -x "FORMAT/PL" \\
        -Ov \\
        -o tmp.${sample_id}.${caller}.vcf \\
        ${vcf_file}

    # 2) Extract header and fix misdeclared INFO fields
    bcftools view -h tmp.${sample_id}.${caller}.vcf | \\
        sed -E 's/(ID=HAPCOMP),Number=A/\\1,Number=1/' | \\
        sed -E 's/(ID=HAPDOM),Number=A/\\1,Number=1/' | \\
        sed -E 's/(ID=X_HIL),Number=A/\\1,Number=1/' | \\
        sed -E 's/(ID=X_IL),Number=A/\\1,Number=1/' \\
        > header_fixed.hdr

    # 3) Reheader with the fixed header
    bcftools reheader \\
        -h header_fixed.hdr \\
        tmp.${sample_id}.${caller}.vcf \\
        > ${sample_id}.${caller}.vcf
    """

    stub:
    """
    touch "${sample_id}.${caller}.vcf"
    """
}
