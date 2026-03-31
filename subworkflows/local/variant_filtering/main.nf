/*
 * Subworkflow: variant_filtering
 * Purpose: Subset variants, extract VCF annotations, plot distributions, and apply final filtering.
 * Inputs:
 *   - ch_vcf: Channel emitting tuple(meta, vcf_file) [bgzipped and indexed]
 *   - ch_fasta: Channel emitting tuple(meta, fasta_file)
 *   - ch_output_dir: Channel emitting output directory path
 *   - ch_annotations: Channel emitting list of annotation format strings (e.g. ["%SOR", "%FS", ...])
 * Outputs:
 *   - snp_vcf: tuple(meta, snp_vcf)
 *   - indel_vcf: tuple(meta, indel_vcf)
 *   - filtered_vcf: tuple(meta, filtered_vcf)
 *   - plot_png: tuple(meta, annotation, plot PNG)
 *   - annotation_txt: tuple(meta, annotation, annotation values text file)
 *   - versions: version info from all modules
 */

include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS__SNP } from '../../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS__INDEL } from '../../../modules/nf-core/gatk4/selectvariants/main.nf'
include { GATK4_VARIANTSTOTABLE } from '../../../modules/nf-core/gatk4/variantstotable/main'

include { UTILS_FINALFILTER      } from '../../../modules/local/utils/finalfilter/main.nf'

include { BCFTOOLS_INDEX      } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_QUERY      } from '../../../modules/local/bcftools/query/main.nf'
include { BCFTOOLS_NORM       } from '../../../modules/nf-core/bcftools/norm/main.nf'
include { BCFTOOLS_FILTER     } from '../../../modules/nf-core/bcftools/filter/main.nf'
include { BCFTOOLS_VIEW       } from '../../../modules/nf-core/bcftools/view/main.nf'
include { VCFTOOLS as VCFTOOLS_VIEW       } from '../../../modules/nf-core/vcftools/main.nf'
include { R_PLOT_DENSITY      } from '../../../modules/local/utils/rplot/main.nf'

workflow VARIANT_FILTERING {

    take:
        ch_vcf_tbi         // tuple(meta, vcf_file)
        ch_fasta       // tuple(meta, fasta_file)
        ch_faidx_index // tuple(meta, faidx_index)
        ch_dict_index   // tuple(meta, dict_index)
        ch_intervals_list // path to intervals.list (optional)

    main:

        // Collect version info
        ch_versions = channel.empty()
        ch_intervals_list_val = ch_intervals_list
        ch_fasta_val = ch_fasta
        ch_faidx_index_val = ch_faidx_index
        ch_dict_index_val = ch_dict_index

        // Each tuple in ch_vcf_tbi is processed independently (one run per merged sample VCF/TBI).
    //TODO: 20 instances of this process (per-chr multi-sample)
        GATK4_SELECTVARIANTS__SNP(
            ch_vcf_tbi,
            ch_intervals_list_val,
            ch_fasta_val,
            ch_faidx_index_val,
            ch_dict_index_val,
            'snp'
        )
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS__SNP.out.versions.first())

        ch_snp_vcf = GATK4_SELECTVARIANTS__SNP.out.vcf
        ch_snp_vcf_tbi = GATK4_SELECTVARIANTS__SNP.out.vcf_tbi

    //TODO: 20 instances of this process (per-chr multi-sample)
        GATK4_SELECTVARIANTS__INDEL(
            ch_vcf_tbi,
            ch_intervals_list_val,
            ch_fasta_val,
            ch_faidx_index_val,
            ch_dict_index_val,
            'indel'
        )
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS__INDEL.out.versions.first())

        list_annotation_types = [
            "[ %DP]",
            "[ %GQ]",
            "%AC",
            "%AD",
            "%AF",
            "%AN",
            "%FS",
            "%GQ",
            "%INFO/AF",
            "%INFO/DP",
            "%MQ",
            "%MQRankSum",
            "%QD",
            "%QUAL",
            "%ReadPosRankSum",
            "%DP",
            "%SOR" ]

// NOTE: GATK4_VARIANTSTOTABLE is an alternative for the bcftools query based approach

    //TODO: 20 instances of this process (per-chr multi-sample)
        GATK4_VARIANTSTOTABLE(
            ch_snp_vcf_tbi,
            list_annotation_types,
            ch_fasta_val,
            ch_faidx_index_val,
            ch_dict_index_val
        )

        ch_versions = ch_versions.mix(GATK4_VARIANTSTOTABLE.out.versions.first())

        // 2. Index SNP VCF for downstream querying and filtering
    //TODO: 20 instances of this process (per-chr multi-sample)
        BCFTOOLS_INDEX(ch_snp_vcf)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions.first())
        // ch_snp_vcf_indexed = BCFTOOLS_INDEX.out.tbi

        // 3. Extract annotation fields and plot distributions
    //TODO: 20 instances of this process (per-chr multi-sample)
        BCFTOOLS_QUERY( BCFTOOLS_INDEX.out.vcf_tbi, list_annotation_types)

        ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())
        // ch_annotation_txt = BCFTOOLS_QUERY.out.txt

/*
        // R_PLOT_DENSITY(ch_annotation_txt.map { meta, txt, annotation -> tuple(meta, txt, annotation) })
*/

    //TODO: 20 instances of this process (per-chr multi-sample)
        UTILS_FINALFILTER(
            ch_snp_vcf_tbi,
            ch_fasta_val,
            ch_faidx_index_val
        )


    emit:
        // snp_vcf        = ch_snp_vcf
        // indel_vcf      = ch_indel_vcf
        // filtered_vcf   = ch_final_filtered_vcf
        // plot_png       = R_PLOT_DENSITY.out.plot_png
        // annotation_txt = BCFTOOLS_QUERY.out.txt
        gvcf_tbi     = UTILS_FINALFILTER.out.gvcf_tbi
        versions       = ch_versions
}
