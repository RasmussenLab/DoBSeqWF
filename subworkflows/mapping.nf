/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAPPING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                    } from '../modules/fastqc'
include { MOSDEPTH                  } from '../modules/mosdepth'
include { MOSDEPTH as RAW_DEPTH     } from '../modules/mosdepth'
include { FLAGSTAT                  } from '../modules/flagstat'
include { ALIGNMENT                 } from '../modules/alignment'
include { MARKDUPLICATES_FAST       } from '../modules/markduplicates_fast'
include { MARKDUPLICATES            } from '../modules/markduplicates'
include { CLEAN                     } from '../modules/clean'
include { ADDREADGROUP              } from '../modules/addreadgroup'
include { INDELQUAL                 } from '../modules/indelqual'
include { INDEX                     } from '../modules/index'
include { INDEX as RAW_INDEX        } from '../modules/index'
include { CRAM                      } from '../modules/cram'
include { CRAMTABLE                 } from '../modules/cramtable'
include { VALIDATE                  } from '../modules/validate'

// Full QC modules
include { HS_METRICS                } from '../modules/hs_metrics'
include { ALIGNMENT_METRICS         } from '../modules/alignment_metrics'
include { GC_METRICS                } from '../modules/gc_metrics'
include { INSERT_SIZE_METRICS       } from '../modules/insert_size_metrics'
include { DUPLICATE_METRICS         } from '../modules/duplicate_metrics'

include { MULTIQC                   } from '../modules/multiqc'

workflow MAPPING {
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

    if (params.doFastqc || params.fullQC) {
        FASTQC(read_ch)
        qc_ch = FASTQC.out.fastqc_zip
    }

    ALIGNMENT(read_ch, reference_genome)

    RAW_INDEX(ALIGNMENT.out.raw_bam_file)
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index, bedfile, "")
    ADDREADGROUP(ALIGNMENT.out.raw_bam_file)

    if (params.fast_markdup) {
        MARKDUPLICATES_FAST(ADDREADGROUP.out.bam_file, "")
        markdup_ch = MARKDUPLICATES_FAST.out.marked_bam_file
        markdup_metrics_ch = MARKDUPLICATES_FAST.out.metrics_file
    } else {
        MARKDUPLICATES(ADDREADGROUP.out.bam_file, false, "")
        markdup_ch = MARKDUPLICATES.out.marked_bam_file
        markdup_metrics_ch = MARKDUPLICATES.out.metrics_file
    }

    if (params.fullQC) {
        HS_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome, bedfile, "")
        ALIGNMENT_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome, "")
        GC_METRICS(ALIGNMENT.out.raw_bam_file, reference_genome, "")
        INSERT_SIZE_METRICS(ALIGNMENT.out.raw_bam_file, "")
        DUPLICATE_METRICS(ALIGNMENT.out.raw_bam_file, "")

        qc_ch = qc_ch.mix(
            RAW_DEPTH.out.region_dist,
            ALIGNMENT_METRICS.out.metrics_file,
            HS_METRICS.out.metrics_file,
            GC_METRICS.out.metrics_file,
            GC_METRICS.out.summary_file,
            INSERT_SIZE_METRICS.out.metrics_file,
            markdup_metrics_ch)
    } else {
        qc_ch = qc_ch.mix(
            RAW_DEPTH.out.region_dist,
            markdup_metrics_ch)
    }
    
    CLEAN(markdup_ch)

    // Add indel quality, ie. BI/BD tags (for LoFreq compatibility)
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome)

    CRAM(INDELQUAL.out.bam_file, reference_genome)
    CRAMTABLE(CRAM.out.cram_info.collect())

    INDEX(INDELQUAL.out.bam_file)
    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile, "deduplicated")

    if (!params.testing) {
        VALIDATE(INDELQUAL.out.bam_file)
    }

    MULTIQC(qc_ch.mix(
        MOSDEPTH.out.region_dist).collect())

    emit:
    bam_file_w_index = INDEX.out.bam_file_w_index
}