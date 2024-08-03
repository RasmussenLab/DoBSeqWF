#!/usr/bin/env nextflow

// mads 2024-02-16

// Use newest nextflow dsl
nextflow.enable.dsl = 2

log.info """\
    ===================================
             D o B S e q - W F
    ===================================
    """
    .stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAPPING                   } from './subworkflows/mapping'
include { MAPPING_UMI               } from './subworkflows/mapping_umi'

include { CALLING                   } from './subworkflows/calling'
include { CALL_TRUTH                } from './subworkflows/call_wgs_truthset'

include { ANNOTATION                } from './subworkflows/annotation'

// Cram conversion
include { BAM                       } from './modules/bam'
include { INDEX                     } from './modules/index'

// Pinpoint methods
include { PILOT_PINPOINT            } from './modules/pilot_pinpoint'
include { PINPOINT                  } from './subworkflows/pinpoint'

// Test
include { TEST                      } from './modules/test'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

pooltable_ch = Channel
    .fromPath(params.pooltable)
    .splitCsv(sep: '\t')
    .map { row -> tuple(row[0], [file(row[1]), file(row[2])]) }

if (params.step == 'calling') {
    cramtable_ch = Channel
        .fromPath(params.cramtable)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }
}

if (params.step == 'pinpoint') {
    vcftable_ch = Channel
        .fromPath(params.vcftable)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2]) }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.reference_genome + "*", checkIfExists: true).collect()
bedfile_ch = Channel.fromPath(params.bedfile, checkIfExists: true).collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.step == 'mapping' || params.step == 'all' || params.step == '') {
        if (params.umi) {
            bam_file_w_index_ch = MAPPING_UMI(pooltable_ch, reference_genome_ch, bedfile_ch)
        } else {
            bam_file_w_index_ch = MAPPING(pooltable_ch, reference_genome_ch, bedfile_ch)
        }
    } else if (params.step == 'calling') {
        BAM(cramtable_ch, reference_genome_ch)
        bam_file_w_index_ch = INDEX(BAM.out.bam_file)
    }

    if (params.step == 'calling' || params.step == 'all' || params.step == '') {
        CALLING(bam_file_w_index_ch, reference_genome_ch, bedfile_ch)
        lofreq_ch = CALLING.out.vcf_lofreq
        gatk_ch = CALLING.out.vcf_hc
    } else if (params.step == 'pinpoint') {
        vcftable_ch
            .branch { sample, vcf_file, caller ->
                GATK: caller == 'GATK'
                lofreq: caller == 'lofreq'
                other: true
            }
            .set { vcftable_ch }
        
        gatk_ch = vcftable_ch.GATK
            .map { sample, vcf_file, caller -> tuple(sample, vcf_file) }
        lofreq_ch = vcftable_ch.lofreq
            .map { sample, vcf_file, caller -> tuple(sample, vcf_file) }
    }

    if (params.step == 'pinpoint' || params.step == 'all' || params.step == '') {
        if (params.pinpoint_method == 'pilot') {
            vcf_ch = gatk_ch
                .mix(lofreq_ch)
                .map { pool_id, vcf_file -> vcf_file }
            
            PILOT_PINPOINT(vcf_ch.collect(), file(params.pooltable), file(params.decodetable))
            pin_ch = PILOT_PINPOINT.out.pinned_variants
        } else if (params.pinpoint_method == 'new') {
            PINPOINT(gatk_ch, file(params.decodetable), reference_genome_ch)
            pin_ch = PINPOINT.out.pinned_variants
        }
    }

    if (params.testing && params.step != 'mapping' && params.step != 'calling' && params.step != 'pinpoint') {
        TEST(pin_ch, file(params.snv_list), params.pinpoint_method)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALTERNATIVE WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Call variants in WGS data from CRAM files
// Requires a tsv with CRAM files
workflow truthset {
    cramtable_ch = Channel
        .fromPath(params.cramtable)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1])) }
    
    CALL_TRUTH(cramtable_ch, reference_genome_ch, bedfile_ch)
}

workflow annotate {
    Channel
        .fromPath(params.vcftable)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), row[2]) }
            .branch { sample, vcf_file, caller ->
                GATK: caller == 'GATK'
                lofreq: caller == 'lofreq'
                other: true
            }
            .set { vcftable_ch }
        
    gatk_ch = vcftable_ch.GATK
        .map { sample, vcf_file, caller -> tuple(sample, vcf_file) }
    lofreq_ch = vcftable_ch.lofreq
        .map { sample, vcf_file, caller -> tuple(sample, vcf_file) }

    ANNOTATION(gatk_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDoBSeq-WF is done!\n" : "Oops .. something went wrong" )
}
