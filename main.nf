#!/usr/bin/env nextflow

// mads 2024-02-16

// Yaml parser
import org.yaml.snakeyaml.Yaml

// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Help text for input and validate parameters
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

log.info """\
    ===================================
             D o B S e q - W F
    ===================================
    """
    .stripIndent()

if (params.help) {
   log.info paramsHelp("nextflow run main.nf")
   exit 0
}
// validateParameters()
log.info paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mapping
include { FASTQC                    } from './modules/fastqc'
include { INDEX as INDEX_RAW        } from './modules/index'
include { MOSDEPTH                  } from './modules/mosdepth'
include { FLAGSTAT                  } from './modules/flagstat'
include { ALIGNMENT                 } from './modules/alignment'
include { MARKDUPLICATES            } from './modules/markduplicates'
include { CLEAN                     } from './modules/clean'
include { ADDREADGROUP              } from './modules/addreadgroup'
include { INDEX                     } from './modules/index'
include { VALIDATE                  } from './modules/validate'

// Variant calling
include { HAPLOTYPECALLER           } from './modules/haplotypecaller'
include { INDELQUAL                 } from './modules/indelqual'
include { LOFREQ                    } from './modules/lofreq'
include { FILTER                    } from './modules/filter'

// Variant pinning
include { PINNING                   } from './modules/pinning'

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFAULT OUTPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

fastqc_ch = Channel.empty()
mosdepth_ch = Channel.empty()
flagstat_ch = Channel.empty()
lofreq_ch = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.reference_genome + "*", checkIfExists: true).collect()
bedfile_ch = Channel.fromPath(params.bedfile).collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAPPING WORKFLOW
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

    if (params.doMosdepth) {
        INDEX_RAW(ALIGNMENT.out.raw_bam_file)
        MOSDEPTH(INDEX_RAW.out.bam_file_w_index, bedfile_ch)
        mosdepth_ch = MOSDEPTH.out.region_dist
    }

    if (params.doFlagstat) {
        FLAGSTAT(ALIGNMENT.out.raw_bam_file)
        flagstat_ch = FLAGSTAT.out.flagstat
    }

    MARKDUPLICATES(ALIGNMENT.out.raw_bam_file)

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    ADDREADGROUP(CLEAN.out.clean_bam_file)

    VALIDATE(ADDREADGROUP.out.bam_file)

    MULTIQC(fastqc_ch.mix(mosdepth_ch, flagstat_ch, MARKDUPLICATES.out.metrics_file).collect())

    emit:
    bam_file = ADDREADGROUP.out.bam_file
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALLING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow calling {
    take:
    bam_file

    main:

    // LoFreq indelqual
    INDELQUAL(bam_file, reference_genome_ch)
    INDEX(INDELQUAL.out.iq_bam_file)

    // LoFreq
    LOFREQ(INDEX.out.bam_file_w_index, reference_genome_ch, bedfile_ch)
    lofreq_ch = LOFREQ.out.vcf_file

    if (params.minAltSupport != 0) {
        FILTER(LOFREQ.out.vcf_file, params.minAltSupport)
        lofreq_ch = FILTER.out.filtered_vcf_file
    }

    // GATK
    HAPLOTYPECALLER(INDEX.out.bam_file_w_index, reference_genome_ch, bedfile_ch)

    // Combine VCF file channels
    HAPLOTYPECALLER.out.vcf_file
        .mix(lofreq_ch)
        .map { pool_id, vcf_file -> vcf_file }
        .set { vcf_file_ch }
    
    // Pin variants
    PINNING(vcf_file_ch.collect(), file(params.pooltable), file(params.decodetable))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    mapping(pooltable_ch)
    calling(mapping.out.bam_file)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDoBSeq-WF is done!\n" : "Oops .. something went wrong" )
}
