/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL TRUTHSET SUBSET SUB-WORKFLOW
        Subsamples alignment to specific region(s) and calls variants.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUBSET                    } from '../modules/subset'
include { ADDREADGROUP              } from '../modules/addreadgroup'
include { MARKDUPLICATES_FAST       } from '../modules/markduplicates_fast'
include { CLEAN                     } from '../modules/clean'
include { INDELQUAL                 } from '../modules/indelqual'
include { INDEX                     } from '../modules/index'
include { BQSR                      } from '../modules/bqsr'
include { APPLY_BQSR                } from '../modules/apply_bqsr'
include { CRAM                      } from '../modules/cram'
include { CRAMTABLE                 } from '../modules/cramtable'
include { MOSDEPTH                  } from '../modules/mosdepth'
include { INDEX_VCF                 } from '../modules/index_vcf'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'

include { HC_TRUTH                  } from '../modules/haplotypecaller_truth'

include { HC_TRUTH_JOINT            } from '../modules/haplotypecaller_truth_joint'
include { GENOMICSDB                } from '../modules/genomicsdb'
include { GENOTYPEGVCF              } from '../modules/genotypegvcf'

include { LOFREQ                    } from '../modules/lofreq'
include { DEEPVARIANT               } from '../modules/deepvariant'

include { MULTIQC                   } from '../modules/multiqc'


workflow CALL_TRUTH {
    take:
    cramtable
    reference_genome
    bedfile

    main:

    // Target regions extraction from truth wgs bam file
    bedfile_bam_extraction_ch = Channel.fromPath(params.bedfile_bam_extraction).collect()

    mills_ch = Channel.fromPath(params.mills + "*", checkIfExists: true).collect()
    g1000_ch = Channel.fromPath(params.g1000 + "*", checkIfExists: true).collect()

    // Subset alignment
    SUBSET(cramtable, reference_genome, bedfile_bam_extraction_ch)

	ADDREADGROUP(SUBSET.out.bam_file)
    MARKDUPLICATES_FAST(ADDREADGROUP.out.bam_file, "")
    CLEAN(MARKDUPLICATES_FAST.out.marked_bam_file)

    // Add indel quality, ie. BI/BD tags
    INDELQUAL(CLEAN.out.clean_bam_file, reference_genome)

    // Base recalibration
    BQSR(INDELQUAL.out.bam_file, reference_genome, mills_ch, g1000_ch)
    
    // Apply recalibration
    APPLY_BQSR(BQSR.out.bqsr_file, reference_genome)

    // Store CRAM files
    CRAM(APPLY_BQSR.out.corrected_bam_file, reference_genome)
    CRAMTABLE(CRAM.out.cram_info.collect())

    // Index
    INDEX(APPLY_BQSR.out.corrected_bam_file)

    // Call variants DeepVariant
    if (params.deepvariant) {
        DEEPVARIANT(INDEX.out.bam_file_w_index, reference_genome, bedfile)
    }

    // Call variants LoFreq
    if (params.lofreq) {
        LOFREQ(INDEX.out.bam_file_w_index, reference_genome, bedfile)
    }

    // Call variants GATK
    if (params.gatk_joint_calling) {
        HC_TRUTH_JOINT(INDEX.out.bam_file_w_index, reference_genome, bedfile)
        GENOMICSDB(HC_TRUTH_JOINT.out.gvcf_file.collect(),
                HC_TRUTH_JOINT.out.gvcf_index.collect(),
                bedfile)
        
        GENOTYPEGVCF(GENOMICSDB.out.gendb, reference_genome)
    } else {
        HC_TRUTH(INDEX.out.bam_file_w_index, reference_genome, bedfile)
        INDEX_VCF(HC_TRUTH.out.vcf_file, 'GATK')
        VARTABLE(INDEX_VCF.out.vcf_w_index, 'GATK')
        NORMALISE_VCF(INDEX_VCF.out.vcf_w_index, reference_genome, 'GATK')
    }

    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile, "deduplicated")

    MULTIQC(MOSDEPTH.out.region_dist.mix(
        MARKDUPLICATES_FAST.out.metrics_file).collect())
}