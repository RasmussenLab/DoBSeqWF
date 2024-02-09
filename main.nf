#!/usr/bin/env nextflow

// mads 2024-02-08

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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

sampletable_ch = Channel
    .fromPath(params.sampletable)
    .splitCsv(sep: '\t')
    .map { row -> tuple(row[0], row[1], row[2]) }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.resourceBase + "/" + params.genomeVersion + ".fna", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow variants {
    
    if (params.doFastqc) {
        FASTQC(sampletable_ch)
    }

    // Add fastp here?

    

}

workflow {
    variants()
}

workflow.onComplete {
    log.info ( workflow.success ? "\n DoBSeq pipeline is done!\n" : "Oops .. something went wrong" )
}
