/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BED_ANNOTATE              } from '../modules/bed_annotate'
include { SNPEFF                    } from '../modules/snpeff'
include { SNPSIFT_CLINVAR           } from '../modules/snpsift_clinvar'
include { SNPSIFT_FILTER            } from '../modules/snpsift_filter'

if (params.bedfile_gene_annotation) {
    bedfile = Channel.fromPath(params.bedfile, checkIfExists: true).collect()
}

snpeff_cache_ch = Channel.fromPath(params.snpeff_cache, checkIfExists: true).collect()
clinvardb_ch = Channel.fromPath(params.clinvar_db + "*", checkIfExists: true).collect()
snpeff_config_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true).collect()

workflow ANNOTATION {
    take:
    vcf_file

    main:
    
    if (params.bedfile_gene_annotation) {
        BED_ANNOTATE(vcf_file, bedfile, 'GATK')
        vcf_ch = BED_ANNOTATE.out.vcf_file
    }

    SNPEFF(vcf_ch, 'GATK', params.snpeff_db, snpeff_cache_ch, snpeff_config_ch)
    SNPSIFT_CLINVAR(SNPEFF.out.snpeff_vcf, 'GATK', clinvardb_ch)
    SNPSIFT_FILTER(SNPSIFT_CLINVAR.out.clinvar_vcf, 'GATK')

    emit:
    annotated_vcf = SNPSIFT_CLINVAR.out.clinvar_vcf
}