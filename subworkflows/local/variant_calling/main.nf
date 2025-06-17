/*
 * Variant calling subworkflow: BAM → per-chromosome BAM → GVCF
 * Inputs:
 *   - ch_bam: tuple(meta, bam file)
 *   - ch_fasta: tuple(meta, fasta file)
 * Outputs:
 *   - gvcf_files: tuple(meta, chr, gvcf file)
 *   - versions: version info from all modules
 */

include { SAMTOOLS_INDEX      } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_IDXSTATS  } from '../../../modules/nf-core/samtools/idxstats/main.nf'
include { SAMTOOLS_VIEW      } from '../../../modules/nf-core/samtools/view/main.nf'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main.nf'

workflow VARIANT_CALLING_PARALLEL {

    take:
        ch_bam_bai      // tuple(meta, bam, bai)
        ch_fasta    // tuple(meta, fasta)

    main:
        // Get chromosome names for each BAM
        SAMTOOLS_IDXSTATS(ch_bam_bai)
        ch_chr_info = SAMTOOLS_IDXSTATS.out.idxstats.map { meta, idxstats_file ->
            def chromosomes = file(idxstats_file).text.readLines()
                .findAll { !it.startsWith('*') }
                .collect { it.split('\t')[0] }
            chromosomes.collect { chr -> tuple(meta, chr) }
        }.flatten()

        // Subset BAMs by chromosome
        ch_bam_chr = ch_chr_info.map { meta, chr ->
            def bam = ch_bam_bai.filter { it[0] == meta }[0][1]
            tuple(meta, bam, chr)
        }
        SAMTOOLS_VIEW(ch_bam_chr)

        // Run HaplotypeCaller per chromosome
        ch_hc_input = SAMTOOLS_VIEW.out.bam.map { meta, chr_bam, chr ->
            tuple(meta, chr_bam, chr, ch_fasta.filter { it[0] == meta }[0][1])
        }
        GATK4_HAPLOTYPECALLER(ch_hc_input)

        // Collect version info
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    emit:
        gvcf_files = GATK4_HAPLOTYPECALLER.out.gvcf
        versions   = ch_versions
}
