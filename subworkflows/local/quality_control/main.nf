#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { FASTQC as FASTQC_PRE } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_POST } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTP }  from '../../../modules/nf-core/fastp/main.nf'

workflow QUALITY_CONTROL {
    take:
        reads_ch

    main:
        // Run FastQC before filtering
        FASTQC_PRE(reads_ch)

        // Run fastp
        FASTP(reads_ch, [], false, false, false)

    emit:
        fastqc_pre_reports_zip = FASTQC_PRE.out.zip
        fastp_json = FASTP.out.json
        fastp_reads = FASTP.out.reads
        versions = FASTQC_PRE.out.versions.mix(FASTP.out.versions)
}
