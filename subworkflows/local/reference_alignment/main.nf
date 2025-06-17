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

include { UNTAR as UNTAR_BWA_INDEX } from '../../../modules/nf-core/untar/main'
include { UNTARFILES           } from '../../../modules/nf-core/untarfiles/main'
include { BWA_INDEX            } from '../../../modules/nf-core/bwa/index/main'

// nf-core subworkflow for reference alignment
include { FASTQ_ALIGN_BWA } from '../../nf-core/fastq_align_bwa/main' 

workflow REFERENCE_ALIGNMENT {
    take:
        ch_fasta  // channel: [ val(meta), path(fasta) ] (optional, if bwa_index is not provided)
        ch_reads    // tuple(meta, fq1, fq2)
        ch_bwa_index // channel: [ val(meta), path(bwa_index) ] (optional, if bwa_index is provided)

    main:

        ch_versions = Channel.empty()

        // Generate from scratch if required

        BWA_INDEX ( ch_fasta )
        ch_bwa_index = BWA_INDEX.out.index
        ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)



        //
        // SUBWORKFLOW: Alignment with BWA & BAM QC
        //
        FASTQ_ALIGN_BWA (
            ch_reads,
            ch_bwa_index,
            false,
            ch_fasta
        )

}

//        // Sort SAM files
//        PICARD_SORTSAM(
//            BWA_MEM.out.bam,
//            'coordinate'
//        )
//        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions.first())
//
//        // Mark and remove duplicates
//        PICARD_MARKDUPLICATES(
//            PICARD_SORTSAM.out.bam,
//            [],  // no reference fasta needed
//            []   // no reference fai needed
//        )
//        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
//
//        // Add or replace read groups
//        ch_readgroup_bam = PICARD_MARKDUPLICATES.out.bam.map { meta, bam ->
//            def new_meta = meta.clone()
//            new_meta.id = meta.id
//            new_meta.sample = meta.id
//            new_meta.library = meta.id
//            new_meta.platform = 'Illumina'
//            new_meta.platform_unit = '1'
//            return [new_meta, bam]
//        }
//
//        PICARD_ADDORREPLACEREADGROUPS(
//            ch_readgroup_bam,
//            [[],[]], // no reference fasta needed
//            [[],[]]  // no reference fai needed
//        )
//        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())
//
//        // Pre-filtering flagstat
//        ch_bam_for_flagstat_pre = PICARD_ADDORREPLACEREADGROUPS.out.bam.map { meta, bam ->
//            [meta, bam, []]  // Add empty index channel
//        }
//        
//        SAMTOOLS_FLAGSTAT_PRE(
//            ch_bam_for_flagstat_pre
//        )
//        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_PRE.out.versions.first())
//
//        // Filter BAM files (properly paired, not unmapped, not secondary, not supplementary)
//        ch_bam_for_view = PICARD_ADDORREPLACEREADGROUPS.out.bam.map { meta, bam ->
//            [meta, bam, []]  // Add empty index channel
//        }
//        
//        SAMTOOLS_VIEW(
//            ch_bam_for_view,
//            [[],[]], // no reference needed
//            [],      // no qname file
//            []       // no index format
//        )
//        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
//
//        // Post-filtering flagstat
//        ch_bam_for_flagstat_post = SAMTOOLS_VIEW.out.bam.map { meta, bam ->
//            [meta, bam, []]  // Add empty index channel
//        }
//        
//        SAMTOOLS_FLAGSTAT_POST(
//            ch_bam_for_flagstat_post
//        )
//        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_POST.out.versions.first())
//
//        // Combine flagstat reports
//        ch_flagstat_combined = SAMTOOLS_FLAGSTAT_PRE.out.flagstat
//            .join(SAMTOOLS_FLAGSTAT_POST.out.flagstat)
//            .map { meta, pre, post -> 
//                [meta, [pre, post]]
//            }
//    emit:
//        bam_files = SAMTOOLS_VIEW.out.bam           // tuple(meta, bam)
//        flagstat_reports = ch_flagstat_combined     // tuple(meta, [pre_flagstat, post_flagstat])
//        metrics = PICARD_MARKDUPLICATES.out.metrics // tuple(meta, metrics)
//        versions = ch_versions                      // versions from all modules
//}