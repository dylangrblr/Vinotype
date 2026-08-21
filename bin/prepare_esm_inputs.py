#!/usr/bin/env python3
"""Prepare per-transcript FASTA and variant input files for ESM scoring.

Reads the filtered annotation TSV produced by SNPSIFT_EXTRACTFIELDS
(columns: CHROM POS REF ALT ANN[*].IMPACT ANN[*].EFFECT ANN[*].GENEID
ANN[*].FEATUREID ANN[*].HGVS_P) and the combined transcript FASTA from
SAMTOOLS_FAIDX_EXTRACT, then emits:

  fasta_dir/{transcript_id}.fasta   — one FASTA per transcript
  var_dir/{transcript_id}.var.txt   — one variant per line: {WT}{POS}{MT}
  {prefix}.mapping.tsv              — original TSV fields + variant_key column
"""
import argparse
import os
import re
import sys
from collections import defaultdict

from Bio import SeqIO


AA3_TO_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*', 'Sec': 'U', 'Pyl': 'O',
}


def parse_hgvs_p(hgvs_p: str) -> tuple | None:
    """Convert p.Ala123Val → ('A', 123, 'V').  Returns None if unparseable."""
    m = re.match(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', hgvs_p.strip())
    if not m:
        return None
    wt = AA3_TO_1.get(m.group(1))
    mt = AA3_TO_1.get(m.group(3))
    if not wt or not mt:
        return None
    return wt, int(m.group(2)), mt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tsv',      required=True, help='Annotation TSV from SNPSIFT_EXTRACTFIELDS')
    parser.add_argument('--fasta',    required=True, help='Combined transcript FASTA from SAMTOOLS_FAIDX_EXTRACT')
    parser.add_argument('--fasta_dir', default='fasta_dir', help='Output dir for per-transcript FASTAs')
    parser.add_argument('--var_dir',   default='var_dir',   help='Output dir for per-transcript variant files')
    parser.add_argument('--prefix',    default='esm',       help='Prefix for mapping TSV')
    args = parser.parse_args()

    os.makedirs(args.fasta_dir, exist_ok=True)
    os.makedirs(args.var_dir,   exist_ok=True)

    seqs = {rec.id: rec for rec in SeqIO.parse(args.fasta, 'fasta')}
    print(f'Transcript sequences loaded: {len(seqs)}', file=sys.stderr)

    # Group variants by transcript; build per-row mapping for later re-annotation
    transcript_variants: dict[str, list[str]] = defaultdict(list)
    mapping_rows: list[dict] = []
    skipped = 0

    with open(args.tsv) as fh:
        raw_header = fh.readline().rstrip('\n')
        col_names  = raw_header.split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            row = dict(zip(col_names, parts))
            transcript_id = row.get('ANN[*].FEATUREID', '').strip()
            hgvs_p        = row.get('ANN[*].HGVS_P',   '').strip()

            parsed = parse_hgvs_p(hgvs_p)
            if not parsed:
                print(f'SKIP unparseable HGVS_P "{hgvs_p}" for {transcript_id}', file=sys.stderr)
                skipped += 1
                continue
            if transcript_id not in seqs:
                print(f'SKIP {transcript_id}: not present in FASTA', file=sys.stderr)
                skipped += 1
                continue

            wt, pos, mt = parsed
            variant_key = f'{wt}{pos}{mt}'
            transcript_variants[transcript_id].append(variant_key)
            mapping_rows.append({**row, 'variant_key': variant_key})

    print(f'Transcripts with variants : {len(transcript_variants)}', file=sys.stderr)
    print(f'Variant-transcript pairs  : {sum(len(v) for v in transcript_variants.values())}', file=sys.stderr)
    if skipped:
        print(f'Skipped                   : {skipped}', file=sys.stderr)

    for tid, variants in transcript_variants.items():
        SeqIO.write(seqs[tid],
                    os.path.join(args.fasta_dir, f'{tid}.fasta'), 'fasta')
        with open(os.path.join(args.var_dir, f'{tid}.var.txt'), 'w') as fh:
            for v in sorted(set(variants)):
                fh.write(v + '\n')

    if mapping_rows:
        mapping_header = list(mapping_rows[0].keys())
        with open(f'{args.prefix}.mapping.tsv', 'w') as fh:
            fh.write('\t'.join(mapping_header) + '\n')
            for row in mapping_rows:
                fh.write('\t'.join(row.get(k, '') for k in mapping_header) + '\n')


if __name__ == '__main__':
    main()
