/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL TRUTHSET SUBSET SUB-WORKFLOW
        Aligns according to GATK and subsets to specific region(s) and calls variants.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGNMENT                 } from '../modules/alignment'
include { BAM                       } from '../modules/bam'
include { SUBSET                    } from '../modules/subset'
include { ADDREADGROUP              } from '../modules/addreadgroup'
include { MARKDUPLICATES_FAST       } from '../modules/markduplicates_fast'
include { CLEAN                     } from '../modules/clean'
include { INDELQUAL                 } from '../modules/indelqual'
include { INDEX                     } from '../modules/index'
include { BQSR                      } from '../modules/bqsr'
include { APPLY_BQSR                } from '../modules/apply_bqsr'
include { CRAM                      } from '../modules/cram'
include { CRAM as CRAM_SUBSET       } from '../modules/cram'
include { CRAMTABLE                 } from '../modules/cramtable'
include { MOSDEPTH                  } from '../modules/mosdepth'

include { HC_TRUTH                  } from '../modules/haplotypecaller_truth'

include { HC_TRUTH_JOINT            } from '../modules/haplotypecaller_truth_joint'
include { GENOMICSDB                } from '../modules/genomicsdb'
include { GENOTYPEGVCF              } from '../modules/genotypegvcf'

include { LOFREQ                    } from '../modules/lofreq'
include { DEEPVARIANT               } from '../modules/deepvariant'

include { INDEX_VCF                 } from '../modules/index_vcf'
include { VARTABLE                  } from '../modules/vartable'
include { NORMALISE_VCF             } from '../modules/normalise_vcf'
include { VCFTABLE                  } from '../modules/vcftable'

include { ANNOTATION                } from '../subworkflows/annotation'

include { MULTIQC                   } from '../modules/multiqc'

duplication_info = Channel.empty()

workflow CALL_TRUTH {
    take:
    reference_genome
    bedfile

    main:

    if (params.wgs_subset_input == "" & !params.subset_wgs) {
        // Target regions extraction from truth wgs bam file
        bedfile_bam_extraction_ch = Channel.fromPath(params.bedfile_bam_extraction).collect()

        fqtable = Channel
            .fromPath(params.fqtable)
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], [file(row[1]), file(row[2])]) }

        mills_ch = Channel.fromPath(params.mills + "*", checkIfExists: true).collect()
        g1000_ch = Channel.fromPath(params.g1000 + "*", checkIfExists: true).collect()

        // Align
        ALIGNMENT(fqtable, reference_genome)

        ADDREADGROUP(ALIGNMENT.out.raw_bam_file)
        MARKDUPLICATES_FAST(ADDREADGROUP.out.bam_file, "")
        CLEAN(MARKDUPLICATES_FAST.out.marked_bam_file)

        // Add indel quality, ie. BI/BD tags
        INDELQUAL(CLEAN.out.clean_bam_file, reference_genome)

        // Base recalibration
        BQSR(INDELQUAL.out.bam_file, reference_genome, mills_ch, g1000_ch)
        
        // Apply recalibration
        APPLY_BQSR(BQSR.out.bqsr_file, reference_genome)

        // Store WGS CRAM files
        CRAM(APPLY_BQSR.out.corrected_bam_file, reference_genome)
        
        // Subset alignment
        SUBSET(APPLY_BQSR.out.corrected_bam_file, reference_genome, bedfile_bam_extraction_ch)
        
        // Store subset CRAM files
        CRAM_SUBSET(SUBSET.out.bam_file, reference_genome)

        // Save CRAMTABLE for easy downstream analysis
        CRAMTABLE(CRAM_SUBSET.out.cram_info.collect())
        bam_ch = SUBSET.out.bam_file
        duplication_info = MARKDUPLICATES_FAST.out.metrics_file.collect()
    } else if (params.subset_wgs) {
        cramtable_ch = Channel
            .fromPath(params.cramtable)
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], file(row[1])) }
        // Target regions extraction from truth wgs bam file
        bedfile_bam_extraction_ch = Channel.fromPath(params.bedfile_bam_extraction).collect()
        SUBSET(cramtable_ch, reference_genome, bedfile_bam_extraction_ch)
        CRAM_SUBSET(SUBSET.out.bam_file, reference_genome)
        CRAMTABLE(CRAM_SUBSET.out.cram_info.collect())
        bam_ch = SUBSET.out.bam_file
    } else {
        cramtable_ch = Channel
            .fromPath(params.wgs_subset_input)
            .splitCsv(sep: '\t')
            .map { row -> tuple(row[0], file(row[1])) }
        BAM(cramtable_ch, reference_genome)
        bam_ch = BAM.out.bam_file
    }

    // Index
    INDEX(bam_ch)

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
        NORMALISE_VCF(INDEX_VCF.out.vcf_w_index, reference_genome, 'GATK')

        if (params.annotate) {
            ANNOTATION(NORMALISE_VCF.out.norm_vcf)
            vcf_ch = ANNOTATION.out.annotated_vcf
        } else {
            vcf_ch = NORMALISE_VCF.out.norm_vcf
        }
        VARTABLE(vcf_ch, 'GATK')
    }

    MOSDEPTH(INDEX.out.bam_file_w_index, bedfile, "deduplicated")

    MULTIQC(MOSDEPTH.out.region_dist.mix(
        duplication_info))
}