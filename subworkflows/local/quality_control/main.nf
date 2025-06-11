#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC as FASTQC_PRE } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_POST } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTP }  from '../../../modules/nf-core/fastp/main.nf'

workflow QUALITY_CONTROL_WF {
    take:
        reads_ch

    main:
        // Run FastQC before filtering
        FASTQC_PRE(reads_ch)

/*
        // Run fastp
        fastp_out = FASTP(reads_ch, [], false, false, false)
        // Run FastQC after filtering
        fastqc_post = FASTQC(fastp_out.out.fastp_outputs.map { tuple(it[0], it[1]) })
*/
    emit:
        fastqc_pre_reports_zip = FASTQC_PRE.out.zip
        versions = FASTQC_PRE.out.versions

/*
        fastp_outputs      = fastp_out.out.fastp_outputs
        fastqc_post_reports = fastqc_post.out.fastqc_reports
        versions = fastqc_pre.out.versions.mix(fastp_out.out.versions, fastqc_post.out.versions)
*/
}
