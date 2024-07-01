/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMI AWARE MAPPING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mapping
include { FASTQC                    } from '../modules/fastqc'
include { MOSDEPTH                  } from '../modules/mosdepth'
include { FLAGSTAT                  } from '../modules/flagstat'
include { ALIGNMENT                 } from '../modules/alignment'
include { MARKDUPLICATES_FAST       } from '../modules/markduplicates_fast'
include { MARKDUPLICATES            } from '../modules/markduplicates'
include { MARKDUPLICATES as RAW_DUP } from '../modules/markduplicates'
include { CLEAN                     } from '../modules/clean'
include { ADDREADGROUP              } from '../modules/addreadgroup'
include { INDELQUAL                 } from '../modules/indelqual'
include { INDEX                     } from '../modules/index'
include { CRAM                      } from '../modules/cram'
include { CRAMTABLE                 } from '../modules/cramtable'
include { VALIDATE                  } from '../modules/validate'

// UMI Mapping
include { DOWNSAMPLE                } from '../modules/downsample'
include { UBAM                      } from '../modules/ubam'
include { ALIGNMENT_UMI             } from '../modules/alignment_umi'
include { EXTRACT_UMI               } from '../modules/extract_umi'
include { FASTQ                     } from '../modules/fastq'
include { MERGEBAM                  } from '../modules/mergebam'
include { HS_METRICS                } from '../modules/hs_metrics'
include { ALIGNMENT_METRICS         } from '../modules/alignment_metrics'
include { GC_METRICS                } from '../modules/gc_metrics'
include { INSERT_SIZE_METRICS       } from '../modules/insert_size_metrics'
include { DUPLICATE_METRICS         } from '../modules/duplicate_metrics'
include { INDEX as RAW_INDEX        } from '../modules/index'
include { FLAGSTAT as RAW_FLAGSTAT  } from '../modules/flagstat'
include { MOSDEPTH as RAW_DEPTH     } from '../modules/mosdepth'
include { GROUP_UMI                 } from '../modules/group_umi'
include { UMI_METRICS               } from '../modules/umi_metrics'
include { CALL_CONSENSUS            } from '../modules/call_consensus'
include { ALIGNMENT_UMI as CONS_ALIGN   } from '../modules/alignment_umi'
include { FASTQ as CONS_FASTQ       } from '../modules/fastq'
include { MERGE_CONSENSUS           } from '../modules/merge_consensus'
include { HS_METRICS as CONS_METRIC } from '../modules/hs_metrics'

include { MULTIQC                   } from '../modules/multiqc'


workflow MAPPING_UMI {
    take:
    pooltable
    reference_genome
    bedfile
    
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
    if (params.doFastqc || params.fullQC) {
        FASTQC(FASTQ.out.fastq_files)
        qc_ch = FASTQC.out.fastqc_zip
    } else {
        qc_ch = Channel.empty()
    }

    // Align UMI extracted FastQ files
    ALIGNMENT_UMI(FASTQ.out.fastq_files, reference_genome)

    // Merge aligned bam with unaligned UMI tagged bam.
    // This is done to keep the UMI tags in the aligned bam file. (removed by bam->fastq conversion) 
    //
    // Join bam channels by sample_id
    ALIGNMENT_UMI.out.raw_bam_file
        .join(EXTRACT_UMI.out.umi_extracted_ubam_file)
        .set { bam_files_joined }
    // Merge files
    MERGEBAM(bam_files_joined, reference_genome)

    RAW_INDEX(MERGEBAM.out.bam_file)
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile, "")
    RAW_DUP(MERGEBAM.out.bam_file, true, "")

    // Group reads by UMI
    GROUP_UMI(MERGEBAM.out.bam_file)

    // Call consensus reads.
    CALL_CONSENSUS(GROUP_UMI.out.grouped_bam_file)
    
    // Convert consensus unaligned bam to FastQ
    CONS_FASTQ(CALL_CONSENSUS.out.consensus_ubam_file)

    // Align consensus
    CONS_ALIGN(CONS_FASTQ.out.fastq_files, reference_genome)

    // Merge consenus ubam and consensus bam
    CONS_ALIGN.out.raw_bam_file
        .join(CALL_CONSENSUS.out.consensus_ubam_file)
        .set { consensus_bam_files_joined }
    MERGE_CONSENSUS(consensus_bam_files_joined, reference_genome)

    // Add read group to final bam file for GATK compatibility.
    ADDREADGROUP(MERGE_CONSENSUS.out.bam_file)

    MARKDUPLICATES(ADDREADGROUP.out.bam_file, true, "deduplicated")

    CLEAN(MARKDUPLICATES.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome)

    CRAM(INDELQUAL.out.bam_file, reference_genome)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)

    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile, "deduplicated")

    VALIDATE(INDELQUAL.out.bam_file)

    if (params.fullQC ) {
        // Metrics before consensus calling
        HS_METRICS(MERGEBAM.out.bam_file, reference_genome, bedfile, "")
        ALIGNMENT_METRICS(MERGEBAM.out.bam_file, reference_genome, "")
        GC_METRICS(MERGEBAM.out.bam_file, reference_genome, "")
        INSERT_SIZE_METRICS(MERGEBAM.out.bam_file, "")
        DUPLICATE_METRICS(MERGEBAM.out.bam_file, "")
        UMI_METRICS(GROUP_UMI.out.grouped_bam_file)
        // After consensus calling (HS metrics)
        CONS_METRIC(ADDREADGROUP.out.bam_file, reference_genome, bedfile, "deduplicated")
        
        qc_ch = qc_ch.mix(
            HS_METRICS.out.metrics_file,
            ALIGNMENT_METRICS.out.metrics_file,
            GC_METRICS.out.metrics_file,
            GC_METRICS.out.summary_file,
            INSERT_SIZE_METRICS.out.metrics_file,
            DUPLICATE_METRICS.out.metrics_file,
            CONS_METRIC.out.metrics_file)
    }
    
    MULTIQC(qc_ch.mix(
        RAW_DEPTH.out.region_dist,
        RAW_DUP.out.metrics_file,
        GROUP_UMI.out.family_metrics,
        MOSDEPTH.out.region_dist,
        MARKDUPLICATES.out.metrics_file).collect())

    emit:
    bam_file_w_index = INDEX.out.bam_file_w_index
}