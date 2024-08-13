/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PINPOINTING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDEX_VCF                 } from '../modules/index_vcf'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { PINPY                     } from '../modules/pinpy'
include { ANNOTATION                } from 'subworkflows/annotation'

workflow PINPOINT {
    take:
    vcf_file
    decode_table
    reference_genome
    
    main:
    INDEX_VCF(vcf_file, 'GATK')
    NORMALISE_VCF(INDEX_VCF.out.vcf_w_index, reference_genome, 'GATK')
    ANNOTATION(NORMALISE_VCF.out.norm_vcf)
    VARTABLE(ANNOTATION.out.annotated_vcf, 'GATK')
    
    ANNOTATION.out.annotated_vcf,
        .map { sample, vcf -> vcf }
        .collect()
        .set { norm_vcfs }
    VARTABLE.out.vartable
        .map { sample, vartable -> vartable }
        .collect()
        .set { vartables }
    
    PINPY(
        norm_vcfs,
        vartables,
        decode_table,
        'GATK')
    emit:
    pinned_variants = PINPY.out.lookup_table
}