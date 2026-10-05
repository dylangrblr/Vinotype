#!/usr/bin/env python3

"""
Analyse missed variants from a VCF file and plot diagnostic figures (AB, MQ, BQSZ)

Usage:
    python3 SI_figure_2_missed_variants.py calls_all.vcf.gz

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cyvcf2 import VCF


def classify_variant(variant):
    alts = variant.ALT
    if not alts or alts == ["."] or alts == ["<*>"]:
        return "AD=0"
    if len(variant.REF) > 1 or any(len(a) > 1 for a in alts):
        return "Complex_variant"
    return "AD>0"


def parse_vcf(vcf_path):
    records = []
    for v in VCF(vcf_path):
        dp4 = v.INFO.get("DP4")
        if dp4 is None:
            continue

        ref_fwd, ref_rev, alt_fwd, alt_rev = dp4
        ref_total = ref_fwd + ref_rev
        alt_total = alt_fwd + alt_rev
        # ref_total == 0 means an undefined ratio (no ref support), not a zero balance
        allele_balance = alt_total / ref_total if ref_total > 0 else np.nan

        records.append(
            {
                "chrom": v.CHROM,
                "pos": v.POS,
                "allele_balance": allele_balance,
                "dp": v.INFO.get("DP"),
                "bqbz": v.INFO.get("BQBZ"),
                "mq": v.INFO.get("MQ"),
                "class": classify_variant(v),
            }
        )

    # a complex/indel locus can appear twice (an ALT record plus a no-ALT placeholder
    # at the same position); collapse these so each locus is counted once
    loci = {}
    for rec in records:
        key = (rec["chrom"], rec["pos"])
        if key not in loci or (loci[key]["class"] == "AD=0" and rec["class"] != "AD=0"):
            loci[key] = rec

    return pd.DataFrame(loci.values())


CLASS_COLORS = {
    "Complex_variant": "#d94332",
    "AD=0": "#0173B2",
    "AD>0": "#EB8942",
}


def plot_ab_scatter(df, ax):
    for cls, group in df.groupby("class"):
        vals = group["allele_balance"].dropna()
        if vals.empty:
            continue
        ax.scatter(
            vals.index,
            vals,
            label=f"{cls} (n={len(vals)})",
            color=CLASS_COLORS.get(cls, "grey"),
            alpha=0.75,
            s=55,
            zorder=3,
        )
    ax.axhline(1, color="black", linestyle="--", linewidth=1.2, label="ALT = REF")
    ax.set_xlabel("Site index")
    ax.set_ylabel("Allele balance (ALT reads / REF reads)")
    ax.set_title("Allele balance per site by variant class (benchmark variants missed by HaplotypeCaller)")
    ax.legend(title="Category")
    ax.grid(axis="y", alpha=0.3)


def _grouped_stacked_hist(df, column, ax, bins):
    # separate hist() calls per class draw from zero and hide each other when bins overlap,
    # so classes are combined into one stacked call sharing the same bin edges.
    data, labels, colors = [], [], []
    for cls in sorted(df["class"].unique()):
        vals = df.loc[df["class"] == cls, column].dropna()
        if vals.empty:
            continue
        data.append(vals)
        labels.append(f"{cls} (n={len(vals)})")
        colors.append(CLASS_COLORS.get(cls, "grey"))
    ax.hist(data, bins=bins, stacked=True, label=labels, color=colors)


def plot_bqbz_hist(df, ax):
    vals = df["bqbz"].dropna()
    bins = np.linspace(vals.min(), vals.max(), 16)
    _grouped_stacked_hist(df, "bqbz", ax, bins)
    ax.axvline(0, color="black", linestyle="--", linewidth=1.2, label="BQBZ = 0")
    ax.set_xlabel("BQBZ")
    ax.set_ylabel("Count")
    ax.set_title("Base quality bias (BQBZ) distribution")
    ax.legend(title="Category")
    ax.grid(axis="y", alpha=0.3)


def plot_mq_hist(df, ax):
    vals = df["mq"].dropna()
    # MQ is integer-valued
    bins = np.arange(vals.min() - 0.5, vals.max() + 1.5, 1)
    _grouped_stacked_hist(df, "mq", ax, bins)
    ax.axvline(40, color="black", linestyle="--", linewidth=1.2, label="MQ = 40")
    ax.set_xlabel("MQ")
    ax.set_ylabel("Count")
    ax.set_title("Mapping quality (MQ) distribution")
    ax.legend(title="Category")
    ax.grid(axis="y", alpha=0.3)


def main(vcf_path):
    df = parse_vcf(vcf_path)
    print(f"Parsed {len(df)} sites\n")
    print(df["class"].value_counts().to_string())

    fig, axes = plt.subplots(3, 1, figsize=(12, 15))
    plot_ab_scatter(df, axes[0])
    plot_bqbz_hist(df, axes[1])
    plot_mq_hist(df, axes[2])

    plt.tight_layout()
    out = "missed_variant_diagnostics.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} calls_all.vcf.gz")
        sys.exit(1)
    main(sys.argv[1])
