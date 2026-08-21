include { ESM_PREPARE_INPUTS } from '../../../modules/local/esm/prepare_inputs/main'
include { ESM_PREDICT        } from '../../../modules/local/esm/main'
include { ESM_ANNOTATE       } from '../../../modules/local/esm/annotate/main'

workflow ESM_SCORING {

    take:
    ch_tsv          // channel: [ val(meta), path(tsv) ]        — from ESM_PREPARATION
    ch_transcripts  // channel: [ val(meta), path(fasta) ]      — from ESM_PREPARATION
    ch_model_cache  // channel: path(model_cache) | value([])   — optional HF cache dir

    main:

    ch_versions = Channel.empty()

    //
    // Prepare protein FASTA and variant files, and emit a mapping TSV
    // that can be mapped back to original genomic coordinates
    //
    ch_prepare_in = ch_tsv.join(ch_transcripts)

    ESM_PREPARE_INPUTS ( ch_prepare_in )
    ch_versions = ch_versions.mix(ESM_PREPARE_INPUTS.out.versions.first())

    //
    // Score every variant with ESM2.  The model is specified via
    // ext.args in modules.config (e.g. --models facebook/esm2_t33_650M_UR50D).
    //
    ESM_PREDICT (
        ESM_PREPARE_INPUTS.out.dirs,
        ch_model_cache
    )
    ch_versions = ch_versions.mix(ESM_PREDICT.out.versions.first())

    //
    // Relink ESM scores to the original CHROM/POS/REF/ALT coordinates.
    //
    ch_annotate_in = ESM_PREPARE_INPUTS.out.mapping.join(ESM_PREDICT.out.scores)

    ESM_ANNOTATE ( ch_annotate_in )
    ch_versions = ch_versions.mix(ESM_ANNOTATE.out.versions.first())

    emit:
    annotated = ESM_ANNOTATE.out.tsv     // channel: [ val(meta), path(tsv) ]
    scores    = ESM_PREDICT.out.scores   // channel: [ val(meta), path(csv) ]
    figures   = ESM_PREDICT.out.figures  // channel: [ val(meta), path(figures/) ]
    versions  = ch_versions              // channel: [ versions.yml ]
}
