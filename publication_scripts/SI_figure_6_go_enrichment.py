#!/usr/bin/env python3
"""Extract HIGH/MODERATE impact gene IDs from a SnpEff-annotated VCF and run gProfiler GO enrichment.

By default, all ANN entries whose impact level is HIGH or MODERATE are
collected. Extracted V5 IDs are written to *_genes_v5.tsv, then
converted to v4.3 IDs (via --conversion-table) and written to *_genes_v4.tsv
before gProfiler enrichment.

Example
-------
python3 SI_figure_6_go_enrichment.py \
    --vcf results/test.ann.vcf.gz \
    --conversion-table V5_to_v4.3_conversions.xlsx \
    --out-prefix test
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from gprofiler import GProfiler

try:
    from goatools.obo_parser import GODag
    _GOATOOLS_AVAILABLE = True
except ImportError:
    _GOATOOLS_AVAILABLE = False


DEFAULT_ORGANISM = "vvinifera"
DEFAULT_SOURCES = ["GO:BP", "GO:MF", "GO:CC"]
DEFAULT_CONVERSION_TABLE: str | None = None
DEFAULT_CONVERSION_SOURCE_COL = "query_id"
DEFAULT_CONVERSION_TARGET_COL = "v4.3"
GENE_RE = re.compile(r"(Vitvi[0-9A-Za-z_.]*?g\d+|LOC\d+)")

SOURCE_COLORS = {
    "GO:BP": "#2C7BB6",
    "GO:MF": "#1A9641",
    "GO:CC": "#FDAE61",
    "KEGG": "#D7191C",
    "Custom GMT": "#6A5ACD",
}


def _get_all_ancestors(go_id: str, godag) -> set[str]:
    ancestors: set[str] = set()
    to_visit = list(godag[go_id].parents) if go_id in godag else []
    while to_visit:
        parent = to_visit.pop()
        if parent.item_id not in ancestors:
            ancestors.add(parent.item_id)
            to_visit.extend(parent.parents)
    return ancestors


def remove_redundant_go_terms(results: pd.DataFrame, obo_path: str | Path) -> pd.DataFrame:
    """Remove ancestor GO terms from enrichment results using goatools DAG traversal."""
    if not _GOATOOLS_AVAILABLE:
        print("Warning: goatools not installed; will skio redundancy filtering.")
        return results

    godag = GODag(str(obo_path), optional_attrs={"relationship"})
    go_terms = results["native"].tolist()
    valid = {t for t in go_terms if t in godag}

    terms_to_remove: set[str] = set()
    for term in valid:
        ancestors = _get_all_ancestors(term, godag)
        terms_to_remove.update(ancestors & valid)

    terms_to_keep = valid - terms_to_remove
    n_removed = len(go_terms) - len(terms_to_keep)
    print(f"Redundancy filter: {len(go_terms)} terms → {len(terms_to_keep)} retained ({n_removed} redundant ancestors removed)")
    return results[results["native"].isin(terms_to_keep)].copy()


def read_id_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Background gene list not found: {path}")

    ids: list[str] = []
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            ids.append(text.split("\t", 1)[0].strip())
    return list(dict.fromkeys(g for g in ids if g))


def load_conversion_map(
    xlsx_path: str | Path,
    source_col: str = DEFAULT_CONVERSION_SOURCE_COL,
    target_col: str = DEFAULT_CONVERSION_TARGET_COL,
    header_row: int = 0,
    skip_rows: list[int] | None = None,
) -> dict[str, str]:
    path = Path(xlsx_path)
    if not path.exists():
        raise FileNotFoundError(f"Conversion table not found: {path}")

    df = pd.read_excel(
        path,
        header=header_row,
        skiprows=skip_rows,
        dtype=str,
        engine="openpyxl",
    )
    df.columns = [str(c).strip() for c in df.columns]

    if source_col in df.columns and target_col in df.columns:
        df = df[[source_col, target_col]].dropna()
    elif source_col in df.columns and len(df.columns) == 1:
        df = df[[source_col]].dropna()
        target_col = source_col
    else:
        missing = {source_col, target_col} - set(df.columns)
        raise KeyError(
            f"Columns not found in conversion table: {missing}. Available: {list(df.columns)}"
        )

    # Build mapping, merging multiple targets for the same source key with semicolons
    mapping: dict[str, str] = {}
    for src, tgt in zip(
        df[source_col].astype(str).str.strip(),
        df[target_col].astype(str).str.strip(),
    ):
        if src in mapping:
            existing = set(mapping[src].split(";"))
            existing.add(tgt)
            mapping[src] = ";".join(sorted(existing))
        else:
            mapping[src] = tgt
    return mapping


def convert_ids(
    input_ids: list[str] | pd.Series,
    conversion_map: dict[str, str] | None,
    keep_unmapped: bool = False,
) -> tuple[list[str], int]:
    if isinstance(input_ids, pd.Series):
        input_ids = input_ids.tolist()

    if not conversion_map:
        unique = list(dict.fromkeys(str(g) for g in input_ids if str(g).strip()))
        return unique, 0

    converted: list[str] = []
    n_unmapped = 0
    for gene_id in input_ids:
        gene_id = str(gene_id).strip()
        if not gene_id:
            continue
        mapped = conversion_map.get(gene_id)
        if mapped:
            mapped_ids = [part.strip() for part in re.split(r"[;,|]", mapped) if part.strip()]
            if mapped_ids:
                converted.extend(mapped_ids)
        elif keep_unmapped:
            converted.append(gene_id)
        else:
            n_unmapped += 1

    return list(dict.fromkeys(converted)), n_unmapped


def split_gene_tokens(value: str) -> list[str]:
    tokens: list[str] = []
    for part in re.split(r"[,&;|]", value):
        token = part.strip().strip("'").strip('"')
        if token and token != ".":
            tokens.append(token)
    return tokens


def first_gene_token(value: str) -> str:
    tokens = split_gene_tokens(value)
    return tokens[0] if tokens else ""


def extract_gene_ids_from_lof(lof_value: object) -> list[str]:
    gene_ids: list[str] = []
    if lof_value is None:
        return gene_ids

    if isinstance(lof_value, (list, tuple)):
        entries = [str(entry).strip() for entry in lof_value if str(entry).strip()]
    else:
        lof_text = str(lof_value).strip()
        if lof_text.startswith("LOF="):
            lof_text = lof_text[4:].strip()
        if lof_text.startswith("(") and lof_text.endswith(")"):
            lof_text = lof_text[1:-1].strip()
        entries = [entry.strip() for entry in lof_text.split(",") if entry.strip()]

    for entry in entries:
        if "|" not in entry:
            token = first_gene_token(entry)
            if token:
                gene_ids.append(token)
            continue

        fields = [part.strip() for part in entry.split("|")]
        if not fields:
            continue

        for candidate in (fields[0], fields[1] if len(fields) > 1 else ""):
            if candidate and GENE_RE.fullmatch(candidate):
                gene_ids.append(candidate)
                break
        else:
            token = first_gene_token(fields[0])
            if token:
                gene_ids.append(token)

    return gene_ids


def collect_gene_ids(
    vcf_path: Path,
    max_variants: int | None,
    include_impacts: list[str] | None = None,
    include_lof: bool = False,
) -> list[str]:
    vcf = VCF(str(vcf_path))
    seen: list[str] = []
    seen_set: set[str] = set()
    impacts_upper = [i.upper() for i in include_impacts] if include_impacts else []

    variant_count = 0
    for variant in vcf:
        variant_count += 1

        # Extract genes from the SnpEff LOF INFO field if requested
        if include_lof:
            lof_value = variant.INFO.get("LOF")
            if lof_value:
                for gene_id in extract_gene_ids_from_lof(lof_value):
                    if gene_id not in seen_set:
                        seen.append(gene_id)
                        seen_set.add(gene_id)

        # Extract genes from SnpEff ANN annotations filtered by impact level.
        if include_impacts:
            ann_value = variant.INFO.get("ANN")
            if ann_value is not None:
                if isinstance(ann_value, (list, tuple)):
                    ann_entries = [str(x) for x in ann_value if str(x).strip()]
                else:
                    ann_entries = [e.strip() for e in str(ann_value).split(",") if e.strip()]

                # Scan ALL annotation entries so no transcript annotation is missed
                for entry in ann_entries:
                    fields = [f.strip() for f in entry.split("|")]
                    # SnpEff ANN format: Allele|Annotation|Impact|Gene_Name|Gene_ID|etc.
                    impact = fields[2] if len(fields) > 2 else ""
                    if not impact or impact.upper() not in impacts_upper:
                        continue

                    # Prefer Gene_ID (index 4), fallback to Gene_Name (index 3)
                    candidate = ""
                    if len(fields) > 4 and fields[4]:
                        candidate = fields[4]
                    elif len(fields) > 3 and fields[3]:
                        candidate = fields[3]
                    else:
                        candidate = first_gene_token(fields[0] if fields else entry)

                    candidate = str(candidate).strip()
                    if candidate and GENE_RE.fullmatch(candidate) and candidate not in seen_set:
                        seen.append(candidate)
                        seen_set.add(candidate)

        if max_variants is not None and variant_count >= max_variants:
            break

    vcf.close()
    return seen


def run_gprofiler(query_genes: list[str], organism: str, background: list[str] | None, sources: list[str]) -> pd.DataFrame:
    if not query_genes:
        return pd.DataFrame()

    gp_client = GProfiler(return_dataframe=True)
    results = gp_client.profile(
        organism=organism,
        query=query_genes,
        background=background,
        sources=sources,
        significance_threshold_method="fdr",
        user_threshold=0.05,
        all_results=False,
        no_evidences=False,
    )
    if results is None:
        return pd.DataFrame()
    return results


def normalise_gprofiler_results(results: pd.DataFrame) -> pd.DataFrame:
    if results.empty:
        return results

    results = results.copy()
    if "name" not in results.columns and "native" in results.columns:
        results["name"] = results["native"]
    if "native" not in results.columns and "name" in results.columns:
        results["native"] = results["name"]
    if "source" not in results.columns:
        results["source"] = "GO:BP"
    if "p_value" in results.columns:
        results["p_value"] = pd.to_numeric(results["p_value"], errors="coerce")
        results["-log10_p_value"] = -np.log10(results["p_value"].clip(lower=1e-300))
    if "intersection_size" in results.columns:
        results["intersection_size"] = pd.to_numeric(results["intersection_size"], errors="coerce").fillna(0)
    elif "intersection" in results.columns:
        overlap = results["intersection"].astype(str).str.extract(r"^(\d+)")[0]
        results["intersection_size"] = pd.to_numeric(overlap, errors="coerce").fillna(0)
    else:
        results["intersection_size"] = 0
    return results


def normalise_bubble_sizes(series: pd.Series, s_min: float = 20.0, s_max: float = 420.0) -> pd.Series:
    vmin, vmax = series.min(), series.max()
    if vmax == vmin:
        return pd.Series([s_min] * len(series), index=series.index)
    return s_min + (series - vmin) / (vmax - vmin) * (s_max - s_min)


def plot_bubble(results: pd.DataFrame, out_prefix: Path, top_n: int = 35) -> None:
    if results.empty:
        return

    plot_df = (
        results.sort_values("-log10_p_value", ascending=False)
        .drop_duplicates(subset="name")
        .head(top_n)
        .copy()
    )
    if plot_df.empty:
        return

    plot_df = plot_df.sort_values("-log10_p_value", ascending=True)
    sizes = normalise_bubble_sizes(plot_df["intersection_size"])

    fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.35)))
    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.scatter(
            row["-log10_p_value"],
            i,
            s=float(sizes.iloc[i]),
            color=SOURCE_COLORS.get(str(row["source"]), "grey"),
            alpha=0.8,
            edgecolors="white",
            linewidth=0.5,
        )

    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels([textwrap.shorten(str(n), width=55) for n in plot_df["name"]], fontsize=8)
    ax.set_xlabel("-log10(p-value)")
    ax.set_title("GO Enrichment Top Terms\n(bubble size proportional to intersection size)")
    ax.axvline(-np.log10(0.05), color="grey", lw=0.8, linestyle="--", alpha=0.7, label="p = 0.05")

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=9, label=s)
        for s, c in SOURCE_COLORS.items()
        if s in results["source"].astype(str).values
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, title="Source", fontsize=8, loc="lower right")

    plt.tight_layout()
    out = f"{out_prefix}_bubble.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def plot_source_breakdown(results: pd.DataFrame, out_prefix: Path) -> None:
    if results.empty or "source" not in results.columns:
        return

    counts = results.groupby("source").size().sort_values(ascending=False)
    if counts.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    colors = [SOURCE_COLORS.get(str(s), "grey") for s in counts.index]
    ax.bar(counts.index.astype(str), counts.values, color=colors, alpha=0.85)
    ax.set_ylabel("Significant Terms")
    ax.set_title("Enriched Terms per Source")
    for i, n in enumerate(counts.values):
        ax.text(i, n + 0.3, str(int(n)), ha="center", fontsize=9)
    plt.tight_layout()
    out = f"{out_prefix}_source_breakdown.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract LOF gene IDs from a VCF and run GO enrichment with gProfiler."
    )
    parser.add_argument("--vcf", required=True, help="Annotated VCF or VCF.GZ file")
    parser.add_argument(
        "--out-prefix",
        default="results/vcf_go_enrichment/vcf_go_enrichment",
        help="Output prefix for TSV files (default: results/vcf_go_enrichment/vcf_go_enrichment)",
    )
    parser.add_argument(
        "--background-genes",
        default=None,
        help="Optional text file with one background gene ID per line",
    )
    parser.add_argument(
        "--conversion-table",
        default=DEFAULT_CONVERSION_TABLE,
        help="Excel workbook mapping V5 ids to v4.3 ids. If omitted, gene IDs are not converted.",
    )
    parser.add_argument(
        "--conversion-source-col",
        default=DEFAULT_CONVERSION_SOURCE_COL,
        help="Source column in the conversion workbook (default: query_id)",
    )
    parser.add_argument(
        "--conversion-target-col",
        default=DEFAULT_CONVERSION_TARGET_COL,
        help="Target column in the conversion workbook (default: v4.3)",
    )
    parser.add_argument(
        "--organism",
        default=DEFAULT_ORGANISM,
        help="gProfiler organism code (default: vvinifera)",
    )
    parser.add_argument(
        "--sources",
        default=",".join(DEFAULT_SOURCES),
        help="Comma-separated gProfiler sources (default: GO:BP,GO:MF,GO:CC)",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=None,
        help="Stop after this many VCF records (useful for quick tests)",
    )
    parser.add_argument(
        "--include-impacts",
        default="HIGH,MODERATE",
        help="Comma-separated SnpEff impact levels to extract from ANN (default: HIGH,MODERATE).",
    )
    parser.add_argument(
        "--include-lof",
        action="store_true",
        default=True,
        help="Also extract genes from the SnpEff LOF INFO field.",
    )
    parser.add_argument(
        "--obo",
        default=None,
        help="Path to go-basic.obo for goatools redundancy filtering of the bubble plot.",
    )
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sources = [source.strip() for source in args.sources.split(",") if source.strip()]

    include_impacts = [s.strip() for s in args.include_impacts.split(",") if s.strip()] if args.include_impacts else None
    gene_ids_v5 = collect_gene_ids(vcf_path, args.max_variants, include_impacts, include_lof=args.include_lof)

    conversion_map: dict[str, str] | None = None
    if args.conversion_table:
        table_path = Path(args.conversion_table)
        if table_path.exists():
            conversion_map = load_conversion_map(
                table_path,
                source_col=args.conversion_source_col,
                target_col=args.conversion_target_col,
            )
        else:
            print(f"Warning: conversion table not found at {table_path}; gene IDs will not be converted.")
    # Write extracted V5 IDs before conversion
    genes_v5_out = out_prefix.parent / f"{out_prefix.name}_genes_v5.tsv"
    pd.DataFrame({"gene_id": gene_ids_v5}).to_csv(genes_v5_out, sep="\t", index=False)

    gene_ids, n_unmapped = convert_ids(gene_ids_v5, conversion_map)
    genes_out = out_prefix.parent / f"{out_prefix.name}_genes_v4.tsv"
    pd.DataFrame({"gene_id": gene_ids}).to_csv(genes_out, sep="\t", index=False)

    background_v5 = read_id_list(Path(args.background_genes)) if args.background_genes else None
    background, bg_unmapped = convert_ids(background_v5, conversion_map) if background_v5 else (None, 0)
    enrichment_df = run_gprofiler(gene_ids, args.organism, background, sources)
    enrichment_df = normalise_gprofiler_results(enrichment_df)
    enrichment_out = out_prefix.parent / f"{out_prefix.name}_gprofiler.tsv"
    enrichment_df.to_csv(enrichment_out, sep="\t", index=False)

    # Bubble plot
    bubble_df = enrichment_df
    if args.obo:
        obo_path = Path(args.obo)
        if obo_path.exists():
            bubble_df = remove_redundant_go_terms(enrichment_df, obo_path)
            nonredundant_out = out_prefix.parent / f"{out_prefix.name}_gprofiler_nonredundant.tsv"
            bubble_df.to_csv(nonredundant_out, sep="\t", index=False)
            print(f"Wrote {len(bubble_df)} non-redundant terms to {nonredundant_out}")
        else:
            print(f"Warning: OBO file not found at {obo_path}; bubble plot will use all enrichment terms.")

    plot_prefix = out_prefix.parent / out_prefix.name
    plot_bubble(bubble_df, plot_prefix)
    plot_source_breakdown(enrichment_df, plot_prefix)

    print(f"Wrote {len(gene_ids_v5)} extracted V5 gene IDs to {genes_v5_out}")
    print(f"Wrote {len(gene_ids)} unique v4.3 gene IDs to {genes_out}")
    print(f"Unmapped V5 query IDs: {n_unmapped}")
    if background_v5:
        print(f"Unmapped V5 background IDs: {bg_unmapped}")
    print(f"Wrote {len(enrichment_df)} enrichment rows to {enrichment_out}")


if __name__ == "__main__":
    main()