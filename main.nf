#!/usr/bin/env nextflow

// mads 2024-02-12

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

include { FASTQC                    } from './modules/fastqc'
include { ALIGNMENT                 } from './modules/alignment'
include { MARKDUPLICATES            } from './modules/markduplicates'
include { CLEAN                     } from './modules/clean'
include { ADDREADGROUP              } from './modules/addreadgroup'
include { INDEX                     } from './modules/index'
include { VALIDATE                  } from './modules/validate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

sampletable_ch = Channel
    .fromPath(params.sampletable)
    .splitCsv(sep: '\t')
    .map { row -> tuple(row[0], [row[1], row[2]]) }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.resourceBase + "/" + params.genomeVersion + ".fna", checkIfExists: true).collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow mapping {
    take:
    sampletable
    
    main:
    if (params.doFastqc) {
        FASTQC(sampletable)
    }

    // Add fastp here?

    ALIGNMENT(sampletable, reference_genome_ch)

    // Add samtools stats and mosdepth

    MARKDUPLICATES(ALIGNMENT.out.raw_bam_file)

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    ADDREADGROUP(CLEAN.out.clean_bam_file)

    INDEX(ADDREADGROUP.out.bam_file)

    VALIDATE(ADDREADGROUP.out.bam_file)
    
    emit:
    bam_file = ADDREADGROUP.out.bam_file
}

workflow {
    mapping(sampletable_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\n DoBSeq pipeline is done!\n" : "Oops .. something went wrong" )
}
