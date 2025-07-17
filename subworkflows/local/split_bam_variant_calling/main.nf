// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

def makeBins(start, end, split_size) {
        def bins = []
        for (int i = start; i <= end; i += split_size) {
            def bin_end = Math.min(i + split_size - 1, end)
            bins << [i, bin_end]
        }
        return bins
    }



workflow SPLIT_BAM_VARIANT_CALLING {

    take:
        ch_bam_bai        // Channel of BAM files
        ch_fasta          // Reference FASTA
        ch_faidx_index    // FASTA index
        ch_dict_index     // Dict index
        ch_intervals_list // Intervals list 

    main:
        ch_versions = Channel.empty()

        ch_bins = Channel.from(makeBins(params.bin_start, params.bin_end, params.bin_split_size))
                        .dump(tag: 'ch_bins', pretty: true)
        

        // 1. Index BAMs if needed
//        SAMTOOLS_INDEX(ch_bam_files)

        // 2. Generate all (bam, bin, chrom) combinations
        ch_bam_bin_chrom = ch_bam_bai
                                .dump(tag: 'ch_bam_bin_chrom', pretty: true)
        //     .combine(ch_bins)
        //     .combine(ch_chroms)
        //     .map { bam, bin, chrom -> tuple(bam, bin, chrom) }

        // 3. Subset BAM by chromosome for each bin
//        SAMTOOLS_VIEW(ch_bam_bin_chrom)

        // 4. Variant calling on each split
//        GATK4_HAPLOTYPECALLER(
//            SAMTOOLS_VIEW.out.bam,
//            ch_fasta,
//            ch_faidx_index,
//            ch_dict_index
//        )

        // 5. Merge VCFs per bin (optional, if needed)
//        GATK4_MERGEVCFS(GATK4_HAPLOTYPECALLER.out.vcf, ch_dict_index)

    emit:
        // merged_vcfs = GATK4_MERGEVCFS.out.vcf
        versions = ch_versions                     // channel: [ versions.yml ]
}
