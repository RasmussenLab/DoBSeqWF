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
include { FILTER_VARIANTS           } from '../modules/filter_variants'
include { ANNOTATION                } from '../subworkflows/annotation'
include { VARTABLE_PINS             } from '../modules/vartable_pins'
include { MERGE_PINS                } from '../modules/mergepins'


if (params.filter) {
    snv_model = Channel.fromPath("$projectDir/assets/filter/logistic_regression_snv_model.joblib", checkIfExists: true).collect()
    indel_model = Channel.fromPath("$projectDir/assets/filter/random_forest_indel_model.joblib", checkIfExists: true).collect()
    model_class = Channel.fromPath("$projectDir/assets/filter/model.py", checkIfExists: true).collect()
}

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
    vcf_ch
        .map { sample, vcf -> vcf }
        .collect()
        .set { norm_vcfs }
    
    model_class.view()
    
    if (params.filter) {
        FILTER_VARIANTS(
            norm_vcfs,
            snv_model,
            indel_model,
            model_class,
            params.filter_method,
            params.filter_indels
        )
    }
    
    VARTABLE(vcf_ch, 'GATK')
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
    MERGE_PINS(PINPY.out.vcf_unique_2d_pins)

    emit:
    pinned_variants = PINPY.out.lookup_table
}