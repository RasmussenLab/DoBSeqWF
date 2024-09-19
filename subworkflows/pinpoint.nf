/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PINPOINTING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDEX_VCF                 } from '../modules/index_vcf'
include { DROP_PL                   } from '../modules/drop_pl'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { PINPY                     } from '../modules/pinpy'
include { ANNOTATION                } from '../subworkflows/annotation'
include { VARTABLE_PINS             } from '../modules/vartable_pins'

workflow PINPOINT {
    take:
    vcf_file
    decode_table
    reference_genome
    
    main:
    INDEX_VCF(vcf_file, 'GATK')
    DROP_PL(INDEX_VCF.out.vcf_w_index, 'GATK')
    NORMALISE_VCF(DROP_PL.out.vcf_file, reference_genome, 'GATK')
    
    if (params.annotate) {
        ANNOTATION(NORMALISE_VCF.out.norm_vcf)
        vcf_ch = ANNOTATION.out.annotated_vcf
    } else {
        vcf_ch = NORMALISE_VCF.out.norm_vcf
    }
    
    VARTABLE(vcf_ch, 'GATK')
    vcf_ch
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
    VARTABLE_PINS(PINPY.out.vcf_unique_pins.flatten())

    emit:
    pinned_variants = PINPY.out.lookup_table
}