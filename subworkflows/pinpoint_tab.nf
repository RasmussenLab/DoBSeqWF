/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PINPOINTING WITH TABULAR OUTPUT SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

import static Utils.indexedFileChannel

include { INDEX_VCF                 } from '../modules/index_vcf'
include { DROP_PL                   } from '../modules/drop_pl'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { BGZIP                     } from '../modules/bgzip'
include { PIN_BASIC                 } from '../modules/pin_basic'
include { UNIQUE_VCF                } from '../modules/unique_vcf'
include { VEP                       } from '../modules/vep'

alphamissense_ch = indexedFileChannel(params.alphamissense_tsv, '.tbi')
clinvar_ch = indexedFileChannel(params.clinvar_db, '.tbi')
danmac_db_ch = indexedFileChannel(params.danmac_db, '.tbi')
blacklist_bed_ch = indexedFileChannel(params.blacklist_bed, '.tbi')
repeatmasker_bed_ch = indexedFileChannel(params.repeatmasker_bed, '.tbi')
gnomad_vcf_ch = indexedFileChannel(params.gnomad_vcf, '.tbi')

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
            params.vep_cache                ?: [],
            params.utr_file                 ?: [],
            alphamissense_ch                ?: [],
            clinvar_ch,
            danmac_db_ch,
            blacklist_bed_ch,
            repeatmasker_bed_ch,
            gnomad_vcf_ch,
            params.loftee_gerp_bw           ?: [],
            params.loftee_human_ancestor    ?: [],
            params.loftee_sqlite            ?: []
        )
    }
}