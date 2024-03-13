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

// UMI Mapping
include { UBAM                      } from './modules/ubam'
include { EXTRACT_UMI               } from './modules/extract_umi'
include { FASTQ                     } from './modules/fastq'
include { MERGEBAM                  } from './modules/mergebam'
include { HSMETRICS                 } from './modules/hsmetrics'
include { MARKDUPLICATES as RAW_DUP } from './modules/markduplicates'
include { INDEX as RAW_INDEX        } from './modules/index'
include { FLAGSTAT as RAW_FLAGSTAT  } from './modules/flagstat'
include { MOSDEPTH as RAW_DEPTH     } from './modules/mosdepth'
include { GROUP_UMI                 } from './modules/group_umi'
include { CALL_CONSENSUS            } from './modules/call_consensus'
include { ALIGNMENT as CONS_ALIGN   } from './modules/alignment'
include { FASTQ as CONS_FASTQ       } from './modules/fastq'
include { MERGE_CONSENSUS           } from './modules/merge_consensus'
include { HSMETRICS as CONS_METRIC  } from './modules/hsmetrics'

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
    UMI AWARE MAPPING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow umi_mapping {
    take:
    pooltable
    
    main:

    // Convert FastQ to unaligned bam file
    UBAM(pooltable_ch)

    // Extract UMI from unaligned bam file
    EXTRACT_UMI(UBAM.out.unaligned_bam_file)

    // Convert unaligned bam file to fastq
    FASTQ(EXTRACT_UMI.out.umi_extracted_ubam_file)

    // Run FastQC on UMI extracted FastQ files

    // Single file channel conversion - run FastQC on single files.
    FASTQ.out.fastq_files
        .flatMap { pool_id, reads ->
            return [tuple(pool_id, reads[0], 1), tuple(pool_id, reads[1], 2)]}
        .set { sep_read_ch }
    FASTQC(sep_read_ch)
    fastqc_ch = FASTQC.out.fastqc_zip

    // Align UMI extracted FastQ files
    ALIGNMENT(FASTQ.out.fastq_files, reference_genome_ch)

    // Merge aligned bam with unaligned UMI tagged bam.
    // This is done to keep the UMI tags in the aligned bam file. (removed by bam->fastq conversion) 
    //
    // Join bam channels by sample_id
    ALIGNMENT.out.raw_bam_file
        .join(EXTRACT_UMI.out.umi_extracted_ubam_file)
        .set { bam_files_joined }
    // Merge files
    MERGEBAM(bam_files_joined, reference_genome_ch)
    
    // Collect Hs metrics - Todo
    // HSMETRICS(MERGEBAM.out.bam_file, reference_genome_ch, bedfile_ch)

    // Metrics before consensus calling
    RAW_DUP(MERGEBAM.out.bam_file, "raw")
    RAW_INDEX(MERGEBAM.out.bam_file)
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile_ch, "raw")
    RAW_FLAGSTAT(MERGEBAM.out.bam_file, "raw")

    // Group reads by UMI
    GROUP_UMI(MERGEBAM.out.bam_file)

    // Call consensus reads.
    CALL_CONSENSUS(GROUP_UMI.out.grouped_bam_file)
    
    // Convert consensus unaligned bam to FastQ
    CONS_FASTQ(CALL_CONSENSUS.out.consensus_ubam_file)

    // Align consensus
    CONS_ALIGN(CONS_FASTQ.out.fastq_files, reference_genome_ch)

    // Merge consenus ubam and consensus bam
    CONS_ALIGN.out.raw_bam_file
        .join(CALL_CONSENSUS.out.consensus_ubam_file)
        .set { consensus_bam_files_joined }
    MERGE_CONSENSUS(consensus_bam_files_joined, reference_genome_ch)

    // Add read group to final bam file for GATK compatibility.
    ADDREADGROUP(MERGE_CONSENSUS.out.bam_file)

    // Repeat HS metrics on final bam file
    // CONS_METRIC(ADDREADGROUP.out.bam_file, reference_genome_ch, bedfile_ch)

    MARKDUPLICATES(ADDREADGROUP.out.bam_file, "")

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    // Metrics after consensus calling
    FLAGSTAT(INDELQUAL.out.bam_file, "")
    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch, "")

    VALIDATE(INDELQUAL.out.bam_file)
    
    MULTIQC(FASTQC.out.fastqc_zip.mix(
        RAW_DEPTH.out.region_dist,
        RAW_DUP.out.metrics_file,
        RAW_FLAGSTAT.out.flagstat,
        MOSDEPTH.out.region_dist,
        FLAGSTAT.out.flagstat,
        MARKDUPLICATES.out.metrics_file).collect())

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
    umi_mapping(pooltable_ch)
}

workflow all {
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
