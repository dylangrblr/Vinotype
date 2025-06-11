/*
 * Alignment subworkflow: QC → Alignment → Post-processing
 * Inputs:
 *   - ch_reads: tuple(meta, fq1, fq2)
 *   - val_index_prefix: BWA index prefix (string or path)
 * Outputs:
 *   - bam_files: tuple(meta, bam file)
 *   - versions: version info from all modules
 */

include { BWA_MEM }    from '../../../modules/nf-core/bwa/mem/main.nf'

workflow REFERENCE_ALIGNMENT {
    take:
        ch_reads

    main:

        BWA_MEM(ch_reads, )
        // PICARD_SORTSAM()
        // PICARD_MARKDUPLICATES()
        // PICARD_ADDORREPLACEREADGROUPS()

        // SAMTOOLS_FLAGSTAT_PRE()
        // SAMTOOLS_VIEW()
        // SAMTOOLS_FLAGSTAT_POST()

    emit:
        versions = BWA_MEM.out.versions
}
