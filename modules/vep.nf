process VEP {
    label 'process_single'
    // Annotate using VEP

    conda "$projectDir/envs/vep/environment.yaml"
    container workflow.containerEngine == 'singularity' ? params.container.singularity.vep : params.container.docker.vep

    publishDir "${params.outputDir}/vep_annotate/", pattern: "annotation*.tsv", mode:'copy'

    input:
    tuple path(vcf_file), path(vcf_index), path(variant_keys)
    path reference_genome
    path cache_dir // optional cache directory
    path utr_annotation // optional UTR annotation file
    tuple path(alphamissense), path(am_index) // optional AlphaMissense file
    tuple path(clinvar), path(clinvar_index) // optional ClinVar VCF
    tuple path(danmac), path(danmar_index) // optional DanMac VCF
    tuple path(blacklist), path(blacklist_index) // optional blacklist BED
    tuple path(repeatmasks), path(repeat_index) // optional repeat masks BED
    tuple path(gnomad), path(gnomad_index) // optional gnomAD VCF
    path loftee_plugin_path, stageAs: "loftee" // optional LOFTEE github repo path
    path loftee_gerp_bigwig // optional LOFTEE GERP BigWig
    path loftee_human_ancestor // optional LOFTEE human ancestor FASTA
    path loftee_conservation // optional LOFTEE conservation SQL

    output:
    path "annotations.tsv"
    path "annotations_w_varid.tsv", emit: annotations_file

    script:
    def db = file(params.reference_genome).getName() + ".fasta"
    def mode = cache_dir ? "--cache --offline --dir_cache ${cache_dir}" : "--database"
    def utr = utr_annotation ? "--plugin UTRAnnotator,file=${utr_annotation}" : ""
    def am = alphamissense ? "--plugin AlphaMissense,file=${alphamissense}" : ""
    def clinv = clinvar ? "--custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" : ""
    def danm = danmac ? "--custom ${danmac},DanMac,vcf,exact,0,AC,AN" : ""
    def blackl = blacklist ? "--custom ${blacklist},blacklist,bed,overlap,0,4" : ""
    def repeatm = repeatmasks ? "--custom ${repeatmasks},repeats,bed,overlap,0,4" : ""
    def gnom = gnomad ? "--custom ${gnomad},gnomAD,vcf,exact,0,AC_joint_nfe,AN_joint_nfe,AF_joint_nfe,AC_genomes_nfe,AN_genomes_nfe,AF_genomes_nfe,AC_exomes_nfe,AN_exomes_nfe,AF_exomes_nfe,grpmax_joint,AC_grpmax_joint,AF_grpmax_joint,grpmax_genomes,AC_grpmax_genomes,AF_grpmax_genomes,grpmax_exomes,AC_grpmax_exomes,AF_grpmax_exomes" : ""
    def loftee_plugin_dir = workflow.containerEngine ? "/plugins" : "./loftee"
    def loft = loftee_gerp_bigwig && loftee_human_ancestor && loftee_conservation ? "--plugin LoF,loftee_path:${loftee_plugin_dir},gerp_bigwig:${loftee_gerp_bigwig},human_ancestor_fa:${loftee_human_ancestor},conservation_file:${loftee_conservation}" : ""
    """
    vep \
        --assembly GRCh38 \
        --format vcf \
        --tab \
        --input_file "${vcf_file}" \
        --output_file "annotations.tsv" \
        --fasta "${db}" \
        --check_ref --dont_skip --show_ref_allele \
        --hgvs --canonical --minimal --mane --biotype --numbers --allele_number --symbol \
        --pick --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \
        --force_overwrite \
        --plugin SpliceRegion,extended \
        ${mode} \
        ${utr} \
        ${am} \
        ${clinv} \
        ${danm} \
        ${blackl} \
        ${repeatm} \
        ${gnom} \
        ${loft}
    
    # Header
    printf 'varid\t' > annotations_w_varid.tsv
    grep '^#[^#]' "annotations.tsv" >> annotations_w_varid.tsv
    # Body
    paste                                       \
        <(grep -F -v '*' "${variant_keys}")    \
        <(grep -v '^#'   "annotations.tsv")             \
        >> annotations_w_varid.tsv
    """

    stub:
    """
    touch "annotations.tsv"
    touch "annotations_w_varid.tsv"
    """
}
