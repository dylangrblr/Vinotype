/*
 * Alignment subworkflow: QC → Alignment → Post-processing
 * Inputs:
 *   - ch_reads: tuple(meta, fq1, fq2)
 *   - val_index_prefix: BWA index prefix (string or path)
 * Outputs:
 *   - bam_files: tuple(meta, bam file)
 *   - versions: version info from all modules
 */

//
// Uncompress and prepare reference genome files
//


include { BWAMEM2_MEM                         } from '../../../modules/nf-core/bwamem2/mem/main'
include { PICARD_SORTSAM                      } from '../../../modules/nf-core/picard/sortsam/main'
include { PICARD_MARKDUPLICATES               } from '../../../modules/nf-core/picard/markduplicates/main'
include { PICARD_ADDORREPLACEREADGROUPS       } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { SAMTOOLS_VIEW      } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FLAGSTAT as  SAMTOOLS_FLAGSTAT_PRE      } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as  SAMTOOLS_FLAGSTAT_POST      } from '../../../modules/nf-core/samtools/flagstat/main'

// SUBWORKFLOWS 
include { INDEXES              } from '../../local/indexes/main'

workflow REFERENCE_ALIGNMENT {
    take:
        ch_fasta  // channel: [ val(meta), path(fasta) ] (optional, if bwa_index is not provided)
        ch_reads    // tuple(meta, fq1, fq2)

    main:

        ch_versions = Channel.empty()

        INDEXES (
            ch_fasta
        )


        BWAMEM2_MEM ( ch_reads, INDEXES.out.bwa_index, ch_fasta, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())


        // Sort SAM files
        PICARD_SORTSAM(
             BWAMEM2_MEM.out.bam,
             'coordinate'
        )
        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions.first())

        // Mark and remove duplicates
        PICARD_MARKDUPLICATES(
            PICARD_SORTSAM.out.bam,
            [],  // no reference fasta needed
            []   // no reference fai needed
        )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())


        PICARD_ADDORREPLACEREADGROUPS(
            PICARD_MARKDUPLICATES.out.bam,
            [[],[]], // no reference fasta needed
            [[],[]]  // no reference fai needed
        )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())

        // Pre-filtering flagstat
        ch_bam_for_flagstat_pre = (PICARD_ADDORREPLACEREADGROUPS.out.bam).join(PICARD_ADDORREPLACEREADGROUPS.out.bai)
        
        SAMTOOLS_FLAGSTAT_PRE(
            ch_bam_for_flagstat_pre
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_PRE.out.versions.first())


        SAMTOOLS_VIEW(
            ch_bam_for_flagstat_pre,
            [[],[]], // no reference needed
            [],      // no qname file
            []       // no index format
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

        // Post-filtering flagstat
        ch_bam_for_flagstat_post = (SAMTOOLS_VIEW.out.bam).join(SAMTOOLS_VIEW.out.bai)

        
        SAMTOOLS_FLAGSTAT_POST(
            ch_bam_for_flagstat_post
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_POST.out.versions.first())

        // Combine flagstat reports
        ch_flagstat_combined = SAMTOOLS_FLAGSTAT_PRE.out.flagstat
            .join(SAMTOOLS_FLAGSTAT_POST.out.flagstat)
            .map { meta, pre, post -> 
                [meta, [pre, post]]
            }
    emit:
        bam_files = SAMTOOLS_VIEW.out.bam           // tuple(meta, bam)
        flagstat_reports = ch_flagstat_combined     // tuple(meta, [pre_flagstat, post_flagstat])
        metrics = PICARD_MARKDUPLICATES.out.metrics // tuple(meta, metrics)
        bwa_index =  INDEXES.out.bwa_index
        dict_index =  INDEXES.out.dict_index
        faidx_index = INDEXES.out.faidx_index
        versions = ch_versions                      // versions from all modules
}