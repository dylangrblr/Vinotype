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

include { GATK4_SELECTVARIANTS } from '../../../modules/nf-core/gatk4/selectvariants/main.nf'
include { BCFTOOLS_INDEX      } from '../../../modules/nf-core/bcftools/index/main.nf'
include { BCFTOOLS_QUERY      } from '../../../modules/nf-core/bcftools/query/main.nf'
include { BCFTOOLS_NORM       } from '../../../modules/nf-core/bcftools/norm/main.nf'
include { BCFTOOLS_FILTER     } from '../../../modules/nf-core/bcftools/filter/main.nf'
include { BCFTOOLS_VIEW       } from '../../../modules/nf-core/bcftools/view/main.nf'
include { VCFTOOLS_VIEW       } from '../../../modules/nf-core/vcftools/view/main.nf'
include { R_PLOT_DENSITY      } from '../../../modules/local/utils/rplot/main.nf'

workflow VARIANT_FILTERING {

    take:
        ch_vcf         // tuple(meta, vcf_file)
        ch_fasta       // tuple(meta, fasta_file)
        ch_output_dir  // output directory
        ch_annotations // list of annotation format strings

    main:
        // 1. Subset SNPs and INDELs using GATK4 SelectVariants
        ch_snp_input = ch_vcf.map { meta, vcf -> tuple(meta, vcf, ch_fasta.first()[1], 'SNP') }
        ch_indel_input = ch_vcf.map { meta, vcf -> tuple(meta, vcf, ch_fasta.first()[1], 'INDEL') }

        GATK4_SELECTVARIANTS(ch_snp_input)
        ch_snp_vcf = GATK4_SELECTVARIANTS.out.vcf

        GATK4_SELECTVARIANTS(ch_indel_input)
        ch_indel_vcf = GATK4_SELECTVARIANTS.out.vcf

        // 2. Index SNP VCF for downstream querying and filtering
        BCFTOOLS_INDEX(ch_snp_vcf)
        ch_snp_vcf_indexed = BCFTOOLS_INDEX.out.tbi

        // 3. Extract annotation fields and plot distributions
        ch_annot_tasks = ch_annotations
            .combine(ch_snp_vcf_indexed)
            .map { annotation, tuple(meta, vcf) -> tuple(meta, vcf, annotation, ch_output_dir) }

        BCFTOOLS_QUERY(ch_annot_tasks.map { meta, vcf, annotation, outdir -> tuple(meta, vcf, annotation) })
        ch_annotation_txt = BCFTOOLS_QUERY.out.txt

        R_PLOT_DENSITY(ch_annotation_txt.map { meta, txt, annotation -> tuple(meta, txt, annotation) })

        // 4. Apply final filtering and normalization
        BCFTOOLS_NORM(ch_snp_vcf_indexed.map { meta, vcf -> tuple(meta, vcf, ch_fasta.first()[1]) })
        ch_norm_vcf = BCFTOOLS_NORM.out.vcf

        BCFTOOLS_FILTER(ch_norm_vcf.map { meta, vcf -> tuple(meta, vcf, "SOR >= 3.5 || FS >= 50.0 || ReadPosRankSum <= -2.0 || MQRankSum <= -3.0 || QD <= 10.0 || COUNT(FORMAT/GQ <= 10) / N_SAMPLES > 0.25 || COUNT(GT[*] == '0/1' & AD[*:1] < 2) / COUNT(GT[*] == '0/1') > 0.0 || COUNT(GT[*] == '1/1' & AD[*:1] < 4) / COUNT(GT[*] == '1/1') > 0.0 || FORMAT/DP < 2 || AC < 2") })
        ch_filtered_vcf = BCFTOOLS_FILTER.out.vcf

        BCFTOOLS_NORM(ch_filtered_vcf.map { meta, vcf -> tuple(meta, vcf, '+both') })
        ch_final_norm_vcf = BCFTOOLS_NORM.out.vcf

        BCFTOOLS_VIEW(ch_final_norm_vcf.map { meta, vcf -> tuple(meta, vcf, '--max-alleles 2') })
        ch_biallelic_vcf = BCFTOOLS_VIEW.out.vcf

        VCFTOOLS_VIEW(ch_biallelic_vcf.map { meta, vcf -> tuple(meta, vcf, '--max-missing 0.8') })
        ch_final_filtered_vcf = VCFTOOLS_VIEW.out.vcf

        // 5. Collect version info
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions)
        ch_versions = ch_versions.mix(R_PLOT_DENSITY.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(VCFTOOLS_VIEW.out.versions)

    emit:
        snp_vcf        = ch_snp_vcf
        indel_vcf      = ch_indel_vcf
        filtered_vcf   = ch_final_filtered_vcf
        plot_png       = R_PLOT_DENSITY.out.plot_png
        annotation_txt = BCFTOOLS_QUERY.out.txt
        versions       = ch_versions
}
