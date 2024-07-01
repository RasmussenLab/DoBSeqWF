/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALLING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INTERVALS                 } from '../modules/intervals'
include { HAPLOTYPECALLER           } from '../modules/haplotypecaller'
include { HAPLOTYPECALLER_JOINT     } from '../modules/haplotypecaller_joint'
include { LOFREQ                    } from '../modules/lofreq'
include { DEEPVARIANT               } from '../modules/deepvariant'
include { CRISP                     } from '../modules/crisp'
include { OCTOPUS                   } from '../modules/octopus'
include { FILTER                    } from '../modules/filter'
include { GENOMICSDB                } from '../modules/genomicsdb'
include { GENOTYPEGVCF              } from '../modules/genotypegvcf'
include { MERGEVCFS                 } from '../modules/mergevcfs'
include { VCFTABLE                  } from '../modules/vcftable'

workflow CALLING {
    take:
    bam_file_w_index
    reference_genome
    bedfile

    main:
    
    if (params.lofreq) {
        LOFREQ(bam_file_w_index, reference_genome, bedfile)
        lofreq_ch = LOFREQ.out.vcf_file
        lofreq_info = LOFREQ.out.vcf_info

        if (params.minAltSupport != 0) {
            FILTER(LOFREQ.out.vcf_file, params.minAltSupport)
            lofreq_ch = FILTER.out.filtered_vcf_file
            lofreq_info = FILTER.out.filtered_vcf_info
        }
    } else {
        lofreq_ch = Channel.empty()
        lofreq_info = Channel.empty()
    }

    if (params.crisp) {
        bam_file_w_index
            .map { pool_id, bam_file, index -> bam_file }
            .set { bam_file_ch }
        CRISP(bam_file_ch.collect(), reference_genome, bedfile)
    }

    if (params.octopus) {
        OCTOPUS(bam_file_w_index, reference_genome, bedfile)
    }

    if (params.runHCParallel) {
        // Split by interval list to multi-thread calling
        Channel
            .fromList(params.intervalList)
            .set { intervalList_ch }
        
        bam_file_w_index
            .combine(intervalList_ch)
            .set { bam_file_w_interval }
        
        INTERVALS(reference_genome, bedfile)
        HAPLOTYPECALLER(bam_file_w_interval, reference_genome, INTERVALS.out.target_list)

        HAPLOTYPECALLER.out.vcf_file
            .groupTuple()
            .map { pool_id, vcf_file -> tuple(pool_id, vcf_file.sort()) }
            .set { vcf_by_intervals }
        
        MERGEVCFS(vcf_by_intervals)
        VCFTABLE(MERGEVCFS.out.vcf_info
            .mix(lofreq_info)
            .collect())
        hc_ch = MERGEVCFS.out.vcf_file
    } else {
        if (params.gatk_joint_calling) {
            HAPLOTYPECALLER_JOINT(bam_file_w_index, reference_genome, bedfile)
            GENOMICSDB(HAPLOTYPECALLER_JOINT.out.gvcf_file.collect(),
                    HAPLOTYPECALLER_JOINT.out.gvcf_index.collect(),
                    bedfile)
            
            GENOTYPEGVCF(GENOMICSDB.out.gendb, reference_genome)
            hc_ch = GENOTYPEGVCF.out.vcf_file
        } else {
            HAPLOTYPECALLER(bam_file_w_index, reference_genome, bedfile)
            VCFTABLE(MERGEVCFS.out.vcf_info
                .mix(lofreq_info)
                .collect())
            hc_ch = HAPLOTYPECALLER.out.vcf_file
        }
    }

    emit:
    vcf_hc = hc_ch
    vcf_lofreq = lofreq_ch
}