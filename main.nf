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
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mapping
include { FASTQC                    } from './modules/fastqc'
include { MOSDEPTH                  } from './modules/mosdepth'
include { FLAGSTAT                  } from './modules/flagstat'
include { ALIGNMENT                 } from './modules/alignment'
include { MARKDUPLICATES            } from './modules/markduplicates'
include { CLEAN                     } from './modules/clean'
include { ADDREADGROUP              } from './modules/addreadgroup'
include { INDELQUAL                 } from './modules/indelqual'
include { INDEX                     } from './modules/index'
include { CRAM                      } from './modules/cram'
include { CRAMTABLE                 } from './modules/cramtable'
include { VALIDATE                  } from './modules/validate'

// Cram conversion
include { BAM                       } from './modules/bam'

// Variant calling
include { HAPLOTYPECALLER           } from './modules/haplotypecaller'
include { LOFREQ                    } from './modules/lofreq'
include { FILTER                    } from './modules/filter'

// Variant pinning
include { PINNING                   } from './modules/pinning'

// Test
include { TEST                      } from './modules/test'

// QC
include { MULTIQC                   } from './modules/multiqc'

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFAULT OUTPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

fastqc_ch = Channel.empty()
mosdepth_ch = Channel.empty()
flagstat_ch = Channel.empty()
lofreq_ch = Channel.empty()
bam_file_ch = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.reference_genome + "*", checkIfExists: true).collect()
bedfile_ch = Channel.fromPath(params.bedfile).collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAPPING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow mapping {
    take:
    pooltable
    
    main:
    if (params.doFastqc) {
        // Single file channel conversion - run FastQC on single files.
        pooltable_ch
            .flatMap { pool_id, reads ->
                return [tuple(pool_id, reads[0], 1), tuple(pool_id, reads[1], 2)]}
            .set { sep_read_ch }
        FASTQC(sep_read_ch)
        fastqc_ch = FASTQC.out.fastqc_zip
    }

    ALIGNMENT(pooltable, reference_genome_ch)

    if (params.doFlagstat) {
        FLAGSTAT(ALIGNMENT.out.raw_bam_file)
        flagstat_ch = FLAGSTAT.out.flagstat
    }

	ADDREADGROUP(ALIGNMENT.out.raw_bam_file)

    MARKDUPLICATES(ADDREADGROUP.out.bam_file)

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    if (params.doMosdepth) {
        MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch)
        mosdepth_ch = MOSDEPTH.out.region_dist
    }

    if (!params.testing) {
        VALIDATE(INDELQUAL.out.bam_file)
    }

    MULTIQC(fastqc_ch.mix(mosdepth_ch, flagstat_ch, MARKDUPLICATES.out.metrics_file).collect())

    emit:
    bam_file_w_index = INDEX.out.bam_file_w_index
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALLING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow calling {
    take:
    bam_file_w_index

    main:
    // LoFreq
    LOFREQ(bam_file_w_index, reference_genome_ch, bedfile_ch)
    lofreq_ch = LOFREQ.out.vcf_file

    if (params.minAltSupport != 0) {
        FILTER(LOFREQ.out.vcf_file, params.minAltSupport)
        lofreq_ch = FILTER.out.filtered_vcf_file
    }

    // GATK
    HAPLOTYPECALLER(bam_file_w_index, reference_genome_ch, bedfile_ch)

    // Combine VCF file channels
    HAPLOTYPECALLER.out.vcf_file
        .mix(lofreq_ch)
        .map { pool_id, vcf_file -> vcf_file }
        .set { vcf_file_ch }
    
    // Pin variants
    PINNING(vcf_file_ch.collect(), file(params.pooltable), file(params.decodetable))

    emit:
    pinned_variants = PINNING.out.pinned_variants
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.step != 'calling') {
        bam_file_w_index_ch = mapping(pooltable_ch)
    } else {
        BAM(cramtable_ch, reference_genome_ch)
        bam_file_w_index_ch = INDEX(BAM.out.bam_file)
    }

    if (params.step != 'mapping') {
        calling(bam_file_w_index_ch)
    }

    if (params.testing && params.step != 'mapping') {
        TEST(calling.out.pinned_variants, file(params.snv_list))
    }
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDoBSeq-WF is done!\n" : "Oops .. something went wrong" )
}
