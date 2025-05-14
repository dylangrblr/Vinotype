/*
 * Alignment subworkflow: QC → Alignment → Post-processing
 * Inputs:
 *   - ch_reads: tuple(meta, fq1, fq2)
 *   - val_index_prefix: BWA index prefix (string or path)
 * Outputs:
 *   - bam_files: tuple(meta, bam file)
 *   - versions: version info from all modules
 */

include { FASTP }      from '../../modules/fastp/main.nf'
include { BWA_MEM }    from '../../modules/bwa_mem/main.nf'
include { SAMTOOLS_SORT } from '../../modules/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from '../../modules/samtools/index/main.nf'

workflow REFERENCE_ALIGNMENT {
    take:
        ch_reads
        val_index_prefix

    main:
        // 1. QC and trimming
        fastp_out = FASTP(ch_reads)
        // 2. Alignment
        bwa_out = BWA_MEM(fastp_out.out.fastp_outputs.map { tuple(it[0], it[1], it[2], val_index_prefix) })
        // 3. Sort BAM
        sort_out = SAMTOOLS_SORT(bwa_out.out.sam.map { tuple(it[0], it[1]) })
        // 4. Index BAM
        index_out = SAMTOOLS_INDEX(sort_out.out.bam.map { tuple(it[0], it[1]) })

        // Collect version info
        versions = fastp_out.out.versions.mix(bwa_out.out.versions).mix(sort_out.out.versions).mix(index_out.out.versions)

    emit:
        bam_files = index_out.out.bam_index.map { tuple(it[0], it[1]) }
        versions = versions
}
