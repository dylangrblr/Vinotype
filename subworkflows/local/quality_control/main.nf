#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_PRE } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_POST } from '../../../modules/nf-core/fastqc/main.nf'
include { FASTP } from '../../../modules/nf-core/fastp/main.nf'

workflow QUALITY_CONTROL {
    take:
    reads_ch

    main:
    // Run FastQC before filtering
    FASTQC_PRE(reads_ch)

    save_trimmed_fail = false
    save_merged = false
    FASTP(
        reads_ch,
        [],
        false,
        save_trimmed_fail,
        save_merged,
    )

    //reports = reports.mix(FASTP.out.json.collect { meta, json -> json })
    //reports = reports.mix(FASTP.out.html.collect { meta, html -> html })

    // Run fastp
    FASTP(reads_ch, [], false, false, false)

    if (params.split_fastq) {
        reads_for_alignment = FASTP.out.reads
            .map { meta, reads ->
                def read_files = reads.sort(false) { a, b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                [meta + [n_fastq: read_files.size()], read_files]
            }
            .transpose()
    }
    else {
        reads_for_alignment = FASTP.out.reads
    }


    emit:
    fastqc_pre_reports_zip = FASTQC_PRE.out.zip
    fastp_json = FASTP.out.json
    fastp_reads = reads_for_alignment
    versions = FASTQC_PRE.out.versions.mix(FASTP.out.versions)
}
