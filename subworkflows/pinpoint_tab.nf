/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PINPOINTING WITH TABULAR OUTPUT SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDEX_VCF                 } from '../modules/index_vcf'
include { DROP_PL                   } from '../modules/drop_pl'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { BGZIP                     } from '../modules/bgzip'
include { PIN_BASIC                 } from '../modules/pin_basic'
include { UNIQUE_VCF                } from '../modules/unique_vcf'
include { VEP                       } from '../modules/vep'

workflow PINPOINT_TAB {
    take:
    vcf_file
    matrix_context
    reference_genome
    
    main:
    INDEX_VCF(vcf_file, 'GATK')
    DROP_PL(INDEX_VCF.out.vcf_w_index, 'GATK')
    NORMALISE_VCF(DROP_PL.out.vcf_file, reference_genome, 'GATK')
    BGZIP(NORMALISE_VCF.out.norm_vcf)
    BGZIP.out.bgzf_file
       .map { sample, vcf, index -> tuple(vcf, index) }
       .collect()
       .set { basic_input }
    PIN_BASIC(basic_input,matrix_context)
    UNIQUE_VCF(basic_input)
    if (params.annotate_vep) {
        VEP(
            UNIQUE_VCF.out.vcf_file,
            reference_genome,
            params.vep_cache           ?: [],  // cache
            params.utr_file            ?: [],  // utr
            params.alphamissense_tsv   ?: [],  // alphamissense
            params.clinvar_db          ?: [],  // clinvar
            params.danmac_db           ?: [],  // danmac
            params.blacklist_bed       ?: [],  // blacklist
            params.repeatmasker_bed    ?: [],  // repeatmasker
            params.gnomad_vcf          ?: [],  // gnomad
            params.loftee_gerp_bw      ?: [],  // loftee gerp
            params.loftee_human_ancestor ?: [],// loftee human ancestor
            params.loftee_sqlite       ?: []   // loftee conservation (sqlite)
        )
    }
}