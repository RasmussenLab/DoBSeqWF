/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PILEUP_CALLING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDEX_VCF                 } from '../modules/index_vcf'
include { DROP_PL                   } from '../modules/drop_pl'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { MATRIX_CONTEXT            } from '../modules/matrix_context'
include { BGZIP                     } from '../modules/bgzip'
include { MPILEUP                   } from '../modules/mpileup'

workflow PILEUP_CALLING {
    take:
    vcf_file
    bam_file_w_index
    pooltable
    decode_table
    reference_genome
    bedfile
    
    main:
    INDEX_VCF(vcf_file, 'GATK')
    DROP_PL(INDEX_VCF.out.vcf_w_index, 'GATK')
    NORMALISE_VCF(DROP_PL.out.vcf_file, reference_genome, 'GATK')
    
    BGZIP(NORMALISE_VCF.out.norm_vcf)
    BGZIP.out.bgzf_file
       .map { sample, vcf, index -> tuple(vcf, index) }
       .collect()
       .set { basic_input }
    MATRIX_CONTEXT(pooltable,[])
    MPILEUP(bam_file_w_index, reference_genome, bedfile)

    

    emit:
    pinned_variants = Channel.empty()
}