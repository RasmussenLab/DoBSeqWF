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
include { MARKDUPLICATES_FAST       } from './modules/markduplicates_fast'
include { MARKDUPLICATES            } from './modules/markduplicates'
include { CLEAN                     } from './modules/clean'
include { ADDREADGROUP              } from './modules/addreadgroup'
include { INDELQUAL                 } from './modules/indelqual'
include { INDEX                     } from './modules/index'
include { CRAM                      } from './modules/cram'
include { CRAMTABLE                 } from './modules/cramtable'
include { VALIDATE                  } from './modules/validate'

// UMI Mapping
include { DOWNSAMPLE                } from './modules/downsample'
include { UBAM                      } from './modules/ubam'
include { ALIGNMENT_UMI             } from './modules/alignment_umi'
include { EXTRACT_UMI               } from './modules/extract_umi'
include { FASTQ                     } from './modules/fastq'
include { MERGEBAM                  } from './modules/mergebam'
include { HS_METRICS                } from './modules/hs_metrics'
include { ALIGNMENT_METRICS         } from './modules/alignment_metrics'
include { GC_METRICS                } from './modules/gc_metrics'
include { INSERT_SIZE_METRICS       } from './modules/insert_size_metrics'
include { DUPLICATE_METRICS         } from './modules/duplicate_metrics'
include { INDEX as RAW_INDEX        } from './modules/index'
include { FLAGSTAT as RAW_FLAGSTAT  } from './modules/flagstat'
include { MOSDEPTH as RAW_DEPTH     } from './modules/mosdepth'
include { GROUP_UMI                 } from './modules/group_umi'
include { UMI_METRICS               } from './modules/umi_metrics'
include { CALL_CONSENSUS            } from './modules/call_consensus'
include { ALIGNMENT_UMI as CONS_ALIGN   } from './modules/alignment_umi'
include { FASTQ as CONS_FASTQ       } from './modules/fastq'
include { MERGE_CONSENSUS           } from './modules/merge_consensus'
include { HS_METRICS as CONS_METRIC } from './modules/hs_metrics'

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

read_ch = Channel.empty()
fastqc_ch = Channel.empty()
mosdepth_ch = Channel.empty()
rawdepth_ch = Channel.empty()
flagstat_ch = Channel.empty()
lofreq_ch = Channel.empty()
bam_file_ch = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE FILE CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

reference_genome_ch = Channel.fromPath(params.reference_genome + "*", checkIfExists: true).collect()

// Target regions bedfile
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

    if (params.downsample > 0) {
        read_ch = DOWNSAMPLE(pooltable, params.downsample)
    } else {
        read_ch = pooltable
    }

    if (params.doFastqc) {
        FASTQC(read_ch)
        fastqc_ch = FASTQC.out.fastqc_zip
    }

    ALIGNMENT(read_ch, reference_genome_ch)

    if (params.doFlagstat) {
        FLAGSTAT(ALIGNMENT.out.raw_bam_file)
        flagstat_ch = FLAGSTAT.out.flagstat
    }

    if (params.doMosdepth) {
        RAW_INDEX(ALIGNMENT.out.raw_bam_file)
        RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile_ch, "raw")
        rawdepth_ch = RAW_DEPTH.out.region_dist
    }

	ADDREADGROUP(ALIGNMENT.out.raw_bam_file)
    MARKDUPLICATES_FAST(ADDREADGROUP.out.bam_file, "raw")
    CLEAN(MARKDUPLICATES_FAST.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    if (params.doMosdepth) {
        MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch, "")
        mosdepth_ch = MOSDEPTH.out.region_dist
    }

    if (!params.testing) {
        VALIDATE(INDELQUAL.out.bam_file)
    }

    MULTIQC(fastqc_ch.mix(
        rawdepth_ch,
        mosdepth_ch,
        MARKDUPLICATES_FAST.out.metrics_file).collect())

    emit:
    bam_file_w_index = INDEX.out.bam_file_w_index
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAPPING SUB-WORKFLOW - FULL QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow mapping_qc {
    take:
    pooltable
    
    main:

    if (params.downsample > 0) {
        read_ch = DOWNSAMPLE(pooltable, params.downsample)
    } else {
        read_ch = pooltable
    }

    FASTQC(read_ch)

    ALIGNMENT(read_ch, reference_genome_ch)
    RAW_INDEX(ALIGNMENT.out.raw_bam_file)
    
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile_ch, "raw")
    HS_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome_ch, bedfile_ch, "raw")
    ALIGNMENT_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome_ch, "raw")
    GC_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome_ch, "raw")
    INSERT_SIZE_METRICS(ALIGNMENT.out.raw_bam_file, "raw")
    DUPLICATE_METRICS(ALIGNMENT.out.raw_bam_file, "raw")
    ADDREADGROUP(ALIGNMENT.out.raw_bam_file)
    MARKDUPLICATES_FAST(ADDREADGROUP.out.bam_file, "raw")
    
    CLEAN(MARKDUPLICATES_FAST.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch, "")

    if (!params.testing) {
        VALIDATE(INDELQUAL.out.bam_file)
    }

    MULTIQC(FASTQC.out.fastqc_zip.mix(
        RAW_DEPTH.out.region_dist,
        ALIGNMENT_METRICS.out.metrics_file,
        HS_METRICS.out.metrics_file,
        GC_METRICS.out.metrics_file,
        GC_METRICS.out.summary_file,
        INSERT_SIZE_METRICS.out.metrics_file,
        MOSDEPTH.out.region_dist,
        MARKDUPLICATES_FAST.out.metrics_file).collect())

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

    if (params.downsample > 0) {
        read_ch = DOWNSAMPLE(pooltable, params.downsample)
    } else {
        read_ch = pooltable
    }

    // Convert FastQ to unaligned bam file
    UBAM(read_ch)

    // Extract UMI from unaligned bam file
    EXTRACT_UMI(UBAM.out.unaligned_bam_file)

    // Convert unaligned bam file to fastq
    FASTQ(EXTRACT_UMI.out.umi_extracted_ubam_file)

    // Run FastQC on UMI extracted FastQ files
    FASTQC(FASTQ.out.fastq_files)

    // Align UMI extracted FastQ files
    ALIGNMENT_UMI(FASTQ.out.fastq_files, reference_genome_ch)

    // Merge aligned bam with unaligned UMI tagged bam.
    // This is done to keep the UMI tags in the aligned bam file. (removed by bam->fastq conversion) 
    //
    // Join bam channels by sample_id
    ALIGNMENT_UMI.out.raw_bam_file
        .join(EXTRACT_UMI.out.umi_extracted_ubam_file)
        .set { bam_files_joined }
    // Merge files
    MERGEBAM(bam_files_joined, reference_genome_ch)

    RAW_INDEX(MERGEBAM.out.bam_file)
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile_ch, "raw")

    MARKDUPLICATES_FAST(MERGEBAM.out.bam_file, "raw")

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

    MARKDUPLICATES(ADDREADGROUP.out.bam_file, "")

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    // Metrics after consensus calling
    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch, "")

    VALIDATE(INDELQUAL.out.bam_file)
    
    MULTIQC(FASTQC.out.fastqc_zip.mix(
        RAW_DEPTH.out.region_dist,
        GROUP_UMI.out.family_metrics,
        MOSDEPTH.out.region_dist,
        MARKDUPLICATES_FAST.out.metrics_file,
        MARKDUPLICATES.out.metrics_file).collect())

    emit:
    bam_file_w_index = INDEX.out.bam_file_w_index
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMI AWARE MAPPING SUB-WORKFLOW - FULL QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow umi_mapping_qc {
    take:
    pooltable
    
    main:

    if (params.downsample > 0) {
        read_ch = DOWNSAMPLE(pooltable, params.downsample)
    } else {
        read_ch = pooltable
    }

    // Convert FastQ to unaligned bam file
    UBAM(read_ch)

    // Extract UMI from unaligned bam file
    EXTRACT_UMI(UBAM.out.unaligned_bam_file)

    // Convert unaligned bam file to fastq
    FASTQ(EXTRACT_UMI.out.umi_extracted_ubam_file)

    // Run FastQC on UMI extracted FastQ files

    // Single file channel conversion - run FastQC on single files.
    FASTQC(FASTQ.out.fastq_files)

    // Align UMI extracted FastQ files
    ALIGNMENT_UMI(FASTQ.out.fastq_files, reference_genome_ch)

    // Merge aligned bam with unaligned UMI tagged bam.
    // This is done to keep the UMI tags in the aligned bam file. (removed by bam->fastq conversion) 
    //
    // Join bam channels by sample_id
    ALIGNMENT_UMI.out.raw_bam_file
        .join(EXTRACT_UMI.out.umi_extracted_ubam_file)
        .set { bam_files_joined }
    // Merge files
    MERGEBAM(bam_files_joined, reference_genome_ch)
    
    // Collect Hs metrics
    HS_METRICS(MERGEBAM.out.bam_file, reference_genome_ch, bedfile_ch, "raw")

    // Metrics before consensus calling
    ALIGNMENT_METRICS(MERGEBAM.out.bam_file, reference_genome_ch, "raw")
    GC_METRICS(MERGEBAM.out.bam_file, reference_genome_ch, "raw")
    INSERT_SIZE_METRICS(MERGEBAM.out.bam_file, "raw")
    DUPLICATE_METRICS(MERGEBAM.out.bam_file, "raw")

    RAW_INDEX(MERGEBAM.out.bam_file)
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile_ch, "raw")

    // Group reads by UMI
    GROUP_UMI(MERGEBAM.out.bam_file)

    // Collect UMI metrics
    UMI_METRICS(GROUP_UMI.out.grouped_bam_file)

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
    CONS_METRIC(ADDREADGROUP.out.bam_file, reference_genome_ch, bedfile_ch, "")

    MARKDUPLICATES(ADDREADGROUP.out.bam_file, "")

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome_ch)

    CRAM(INDELQUAL.out.bam_file, reference_genome_ch)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    // Metrics after consensus calling
    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile_ch, "")

    VALIDATE(INDELQUAL.out.bam_file)
    
    MULTIQC(FASTQC.out.fastqc_zip.mix(
        RAW_DEPTH.out.region_dist,
        ALIGNMENT_METRICS.out.metrics_file,
        HS_METRICS.out.metrics_file,
        GC_METRICS.out.metrics_file,
        GC_METRICS.out.summary_file,
        INSERT_SIZE_METRICS.out.metrics_file,
        GROUP_UMI.out.family_metrics,
        DUPLICATE_METRICS.out.metrics_file,
        CONS_METRIC.out.metrics_file,
        MOSDEPTH.out.region_dist,
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

workflow umi {
    umi_mapping(pooltable_ch)
}

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
