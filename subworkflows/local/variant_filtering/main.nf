include { BCFTOOLS_QUERY } from '../../../modules/nf-core/bcftools/query/main.nf'
include { R_PLOT_DENSITY } from '../../../modules/local/utils/rplot/main.nf'

workflow VARIANT_FILTERING {
    /*
     * Subworkflow for extracting VCF annotation fields and plotting their distributions.
     * Inputs:
     *   - meta:      sample metadata (e.g. [id: 'sample1'])
     *   - vcf_file:  input VCF file (bgzipped and indexed)
     *   - output_dir: directory for outputs
     *   - annotations: list of annotation format strings (e.g. ["%SOR", "%FS", ...])
     * Outputs:
     *   - plot_png: tuple(meta, annotation, plot PNG)
     *   - annotation_txt: tuple(meta, annotation, annotation values text file)
     *   - versions: version info from all modules
     */

    take:
        meta
        vcf_file
        output_dir
        annotations

    main:
        // Create a channel of (meta, vcf_file, annotation, output_dir) for each annotation
        ch_annot = annotations.map { annotation -> tuple(meta, vcf_file, annotation, output_dir) }

        // Extract annotation values
        bcftools_results = BCFTOOLS_QUERY(ch_annot)

        // Plot density for each annotation
        plot_results = R_PLOT_DENSITY(bcftools_results.out.annotation_txt.map { tuple(it[0], it[1], it[2]) })

        // Collect version info
        versions = bcftools_results.out.annotation_txt.map { it[3] }
            .mix(plot_results.out.plot_png.map { it[3] })

    emit:
        plot_png = plot_results.out.plot_png.map { tuple(it[0], it[1], it[2]) }
        annotation_txt = bcftools_results.out.annotation_txt.map { tuple(it[0], it[1], it[2]) }
        versions = versions
}
