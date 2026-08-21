#!/usr/bin/env python3
"""Annotate ESM scores back to the original genomic variant coordinates.

Joins the mapping TSV produced by prepare_esm_inputs.py (which carries the
original VCF fields CHROM, POS, REF, ALT, … plus a variant_key column) with
the scores CSV from esm_vep.py (which carries protein, variant, score, …).

Join key: transcript_id (derived from the 'protein' column by stripping the
.fasta suffix) combined with the variant notation (e.g. A123V).
"""
import argparse
import sys

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mapping', required=True, help='Mapping TSV from prepare_esm_inputs.py')
    parser.add_argument('--scores',  required=True, help='ESM scores CSV from esm_vep.py')
    parser.add_argument('--prefix',  default='esm',  help='Output file prefix')
    args = parser.parse_args()

    mapping = pd.read_csv(args.mapping, sep='\t', dtype=str)
    scores  = pd.read_csv(args.scores)

    # Derive transcript ID from the 'protein' basename by stripping .fasta
    scores['_transcript_id'] = scores['protein'].str.replace(r'\.fasta$', '', regex=True)

    # Composite join key avoids false matches when the same variant appears
    # in multiple transcripts
    mapping['_join_key'] = (
        mapping['ANN[*].FEATUREID'].str.strip() + '|' +
        mapping['variant_key'].str.strip()
    )
    scores['_join_key'] = (
        scores['_transcript_id'].str.strip() + '|' +
        scores['variant'].str.strip()
    )

    score_cols = ['score', 'model', 'interpretation', 'ref_aa', 'alt_aa', 'position']
    score_cols = [c for c in score_cols if c in scores.columns]

    merged = mapping.merge(
        scores[['_join_key'] + score_cols],
        on='_join_key',
        how='left',
    ).drop(columns=['_join_key'])

    n_total  = len(merged)
    n_scored = merged['score'].notna().sum() if 'score' in merged.columns else 0
    print(f'Total variant rows   : {n_total}',              file=sys.stderr)
    print(f'Rows with ESM scores : {n_scored}',             file=sys.stderr)
    print(f'Unmatched rows       : {n_total - n_scored}',   file=sys.stderr)

    out_path = f'{args.prefix}.esm_annotated.tsv'
    merged.to_csv(out_path, sep='\t', index=False)
    print(f'Saved → {out_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
