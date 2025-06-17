#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Variant calling subworkflow: BAM → GVCF → Joint Genotyping
 * Inputs:
 *   - ch_bam: tuple(meta, bam file)
 *   - val_ref_fasta: reference fasta file
 * Outputs:
 *   - vcf_files: tuple(meta, vcf file)
 *   - versions: version info from all modules
 */

include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main.nf'
include { GATK4_GENOTYPEGVCFS }   from '../../../modules/nf-core/gatk4/genotypegvcfs/main.nf'
include { SAMTOOLS_FAIDX }        from '../../../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_DICT }        from '../../../modules/nf-core/samtools/dict/main.nf'

workflow VARIANT_CALLING {
    take:
        ch_bam
        ch_fasta

    main:

        SAMTOOLS_FAIDX(ch_fasta)
        SAMTOOLS_DICT(ch_fasta)

        // 1. HaplotypeCaller per sample
        gvcf_out = GATK4_HAPLOTYPECALLER(ch_bam.map { tuple(it[0], it[1], val_ref_fasta) })

        GATK4_HAPLOTYPECALLER
        PICARD_SPLITSAMBYNUMBEROFREADS
        GATK4_MERGEVCFS


    emit:
        vcf_files = genotype_out.out.vcf_files.map { tuple(it[0], it[1]) }
        versions = versions
}
