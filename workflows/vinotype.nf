/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_vinotype_pipeline'

include { QUALITY_CONTROL } from '../subworkflows/local/quality_control'
include { REFERENCE_ALIGNMENT } from '../subworkflows/local/reference_alignment'
include { VARIANT_CALLING } from '../subworkflows/local/variant_calling'
include { JOINT_GENOTYPING } from '../subworkflows/local/joint_genotyping'
include { VARIANT_FILTERING } from '../subworkflows/local/variant_filtering'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VINOTYPE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta      // channel: fasta file read in from --fasta
    ch_intervals_list // channel: intervals list file read in from --intervals_list

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    ch_samplesheet_fastq = ch_samplesheet
                        .filter { meta, files -> meta.data_type == 'fastq' }
                        .dump(tag: 'ch_samplesheet_fastq', pretty: true)



    QUALITY_CONTROL (
        ch_samplesheet_fastq
    )
    ch_multiqc_files  = ch_multiqc_files.mix(QUALITY_CONTROL.out.fastqc_pre_reports_zip.collect{it[1]})
    ch_multiqc_files  = ch_multiqc_files.mix(QUALITY_CONTROL.out.fastp_json.collect{it[1]})
    ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions.first())

    //
    // SUBWORKFLOW: REFERENCE_ALIGNMENT
    //

    ch_samplesheet_bam = ch_samplesheet
                        .filter { meta, files -> meta.data_type == 'bam_genozip' }
                        .dump(tag: 'ch_samplesheet_bam', pretty: true)


    REFERENCE_ALIGNMENT (
        ch_samplesheet_bam,
        ch_fasta,
        QUALITY_CONTROL.out.fastp_reads
    )
    // ch_multiqc_files = ch_multiqc_files.mix()
    // ch_versions = ch_versions.mix(REFERENCE_ALIGNMENT.out.versions.first())

    //
    // SUBWORKFLOW: VARIANT_CALLING
    //


    ch_samplesheet_gvcf = ch_samplesheet
                        .filter { meta, files -> meta.data_type == 'gvcf_genozip' }
                        .dump(tag: 'ch_samplesheet_gvcf', pretty: true)


    VARIANT_CALLING (
        ch_samplesheet_gvcf,
        REFERENCE_ALIGNMENT.out.bam_bai,
        ch_fasta,
        REFERENCE_ALIGNMENT.out.faidx_index,
        REFERENCE_ALIGNMENT.out.dict_index,
    )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    //
    // SUBWORKFLOW: JOINT_GENOTYPING
    //
    JOINT_GENOTYPING (
        VARIANT_CALLING.out.gvcf_tbi_files,
        ch_intervals_list,
        ch_fasta,
        REFERENCE_ALIGNMENT.out.faidx_index,
        REFERENCE_ALIGNMENT.out.dict_index,
        VARIANT_CALLING.out.gvcf_files,
        VARIANT_CALLING.out.tbi_files,
    )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    //
    // SUBWORKFLOW: JOINT_GENOTYPING
    //
    VARIANT_FILTERING (
        JOINT_GENOTYPING.out.gvcf_tbi_files,
        ch_fasta,
        REFERENCE_ALIGNMENT.out.faidx_index,
        REFERENCE_ALIGNMENT.out.dict_index,
        ch_intervals_list,
    )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())




    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'vinotype_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
