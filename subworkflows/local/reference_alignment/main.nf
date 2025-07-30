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
include { GENOZIP_COMPRESS                    } from '../../../modules/local/genozip/compress/main'
include { GENOZIP_DECOMPRESS                  } from '../../../modules/local/genozip/decompress/main'
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
        ch_samplesheet_bam // channel: samplesheet read in from --input
        ch_reads    // tuple(meta, fq1, fq2)

    main:

        ch_versions = Channel.empty()

        INDEXES (
            ch_fasta
        )


        BWAMEM2_MEM ( ch_reads, INDEXES.out.bwa_index, ch_fasta, false )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

        // ch_test_bam = Channel.fromList([
        //         [[id:"Vs_ITA27_dong", single_end:false], "/home/opc/VINOTYPE/work/3b/eab73ea3a7c182b9b747899b7595a9/Vs_ITA27_dong.bam"],
        //         [[id:"Vv_D25_ORLOVI", single_end:false], "/home/opc/VINOTYPE/work/a4/8069d93e6425a45175b11f9fb1e363/Vv_D25_ORLOVI.bam"]
        //     ])

        // Sort SAM files
        PICARD_SORTSAM(
              //ch_test_bam,
              BWAMEM2_MEM.out.bam,
             'coordinate'
        )
        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions.first())


        PICARD_ADDORREPLACEREADGROUPS(
            PICARD_SORTSAM.out.bam,
        )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions.first())


        // Mark and remove duplicates
        PICARD_MARKDUPLICATES(
            PICARD_ADDORREPLACEREADGROUPS.out.bam
        )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())


        // Pre-filtering flagstat
        ch_bam_for_flagstat_pre = (PICARD_MARKDUPLICATES.out.bam).join(PICARD_MARKDUPLICATES.out.bai)
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

        //NOTE: @abhi18av OPTIONALLY IMPLEMENT THIS POST-FILTERING FLAGSTAT
        ch_bam_for_flagstat_post=SAMTOOLS_VIEW.out.bam_bai
        SAMTOOLS_FLAGSTAT_POST(
            ch_bam_for_flagstat_post
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_POST.out.versions.first())
        // // Combine flagstat reports
        // ch_flagstat_combined = SAMTOOLS_FLAGSTAT_PRE.out.flagstat
        //     .join(SAMTOOLS_FLAGSTAT_POST.out.flagstat)
        //     .map { meta, pre, post ->
        //         [meta, [pre, post]]
        //     }
        // ch_flagstat_combined.dump(tag: 'ch_flagstat_combined')


        //SAMTOOLS_VIEW.out.bam_bai.dump(tag:'SAMTOOLS_VIEW.out.bam_bai')

        if(params.use_genozip && params.genozip_license_file) {

            if(params.genozip_compress_bams ) {
                        GENOZIP_COMPRESS(
                            SAMTOOLS_VIEW.out.bam,
                            params.genozip_license_file,
                            INDEXES.out.fasta_genozip
                        )
                        ch_versions = ch_versions.mix(GENOZIP_COMPRESS.out.versions.first())
                    }


            // If samplesheet has GENOZIP BAM files, decompress them
            GENOZIP_DECOMPRESS(
                ch_samplesheet_bam,
                params.genozip_license_file,
                INDEXES.out.fasta_genozip
            )
            ch_versions = ch_versions.mix(GENOZIP_DECOMPRESS.out.versions.first())

            ch_bam_bai = SAMTOOLS_VIEW.out.bam_bai.mix(GENOZIP_DECOMPRESS.out.bam_bai)
        }


   emit:
        // bam_files = SAMTOOLS_VIEW.out.bam           // tuple(meta, bam)
        // metrics = PICARD_MARKDUPLICATES.out.metrics // tuple(meta, metrics)
        //flagstat_reports = ch_flagstat_combined     // tuple(meta, [pre_flagstat, post_flagstat])
        bam_bai =  ch_bam_bai
        bwa_index =  INDEXES.out.bwa_index
        dict_index =  INDEXES.out.dict_index
        faidx_index = INDEXES.out.faidx_index
        versions = ch_versions                      // versions from all modules
}
