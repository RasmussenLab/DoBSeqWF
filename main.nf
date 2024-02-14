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

// Mapping
include { FASTQC                    } from './modules/fastqc'
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

sampletable_ch = Channel
    .fromPath(params.sampletable)
    .splitCsv(sep: '\t')
    .map { row -> tuple(row[0], [file(row[1]), file(row[2])]) }

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
    sampletable
    
    main:
    if (params.doFastqc) {
        // Single file channel conversion - run FastQC on single files.
        sampletable_ch
            .flatMap { sample_id, reads ->
                return [tuple(sample_id, reads[0], 1), tuple(sample_id, reads[1], 2)]}
            .set { sep_read_ch }
        FASTQC(sep_read_ch)
    }

    // Add fastp here?

    ALIGNMENT(sampletable, reference_genome_ch)

    // Add samtools stats and mosdepth

    MARKDUPLICATES(ALIGNMENT.out.raw_bam_file)

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    ADDREADGROUP(CLEAN.out.clean_bam_file)

    VALIDATE(ADDREADGROUP.out.bam_file)
    
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

    // LoFreq
    if (params.doLofreq) {
        INDELQUAL(bam_file, reference_genome_ch)
        
        bam_file_w_idx = INDEX(INDELQUAL.out.iq_bam_file)
        
        LOFREQ(bam_file_w_idx, reference_genome_ch, bedfile_ch)
    } else {
        bam_file_w_idx = INDEX(bam_file)
    }

    // GATK
    HAPLOTYPECALLER(bam_file_w_idx, reference_genome_ch, bedfile_ch)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    mapping(sampletable_ch)
    calling(mapping.out.bam_file)
}