#!/usr/bin/env python3
"""
Compare the per-variant heterozygosity ratio [HET / (HET + HOM ALT)] between
two VCF datasets and visualise distributions as paired violin plots with a
Mann-Whitney U significance test.

Usage:
    python plot_het_ratio_comparison.py \\
        -a datasetA.vcf.gz -b datasetB.vcf.gz \\
        --label-a DatasetA --label-b DatasetB \\
        -o het_ratio_comparison.png --stats het_ratio_stats.tsv
"""

# -- imports --

from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
import argparse
import warnings
import os

# -- helpers (adapted from figure_6_plot_variant_allele_freq.py) --

def compute_het_ratio(genotypes):
    """HET / (HET + HOM ALT) per variant. Returns None when no alt-carrying samples exist."""
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    het     = sum(gt == 1 for gt in genotypes)
    hom_alt = sum(gt == 3 for gt in genotypes)
    total_alt = het + hom_alt
    return het / total_alt if total_alt > 0 else None


def collect_het_ratios(vcf_path):
    """Return a list of per-variant het ratios for every biallelic site."""
    ratios = []
    for variant in VCF(vcf_path):
        ratio = compute_het_ratio(variant.gt_types)
        if ratio is not None:
            ratios.append(ratio)
    return ratios


# -- statistical test (borrowed from SI_figure_2_variant_quality_metrics_violin.py) --

def mw_test(a, b):
    """Two-tailed Mann-Whitney U; returns (U-statistic, p-value)."""
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = mannwhitneyu(a, b, alternative="two-sided", use_continuity=True)
    return res.statistic, res.pvalue


def sig_label(p):
    if np.isnan(p):
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# -- plot --

COLOR_A = "#4C72B0"
COLOR_B = "#DD8452"


def _stat_box(ax, arr, color, x_anchor):
    """Summary stats inset (borrowed from SI_figure_2_variant_quality_metrics_violin.py)."""
    if len(arr) == 0:
        return
    txt = (
        f"n={len(arr):,}\n"
        f"med={np.median(arr):.3f}\n"
        f"IQR [{np.percentile(arr, 25):.3f}–{np.percentile(arr, 75):.3f}]"
    )
    ax.text(
        x_anchor, 0.97, txt,
        transform=ax.transAxes,
        fontsize=7,
        verticalalignment="top",
        color=color,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor=color, alpha=0.75),
    )


def plot_het_ratio(ratios_a, ratios_b, label_a, label_b, output_path, stats_path):
    arr_a = np.array(ratios_a)
    arr_b = np.array(ratios_b)

    u_stat, p = mw_test(arr_a.tolist(), arr_b.tolist())

    fig, ax = plt.subplots(figsize=(6, 7))

    parts_data, positions, colors = [], [], []
    if len(arr_a):
        parts_data.append(arr_a)
        positions.append(1)
        colors.append(COLOR_A)
    if len(arr_b):
        parts_data.append(arr_b)
        positions.append(2)
        colors.append(COLOR_B)

    if parts_data:
        vp = ax.violinplot(
            parts_data,
            positions=positions,
            widths=0.6,
            showmedians=True,
            showextrema=False,
        )
        for body, color in zip(vp["bodies"], colors):
            body.set_facecolor(color)
            body.set_alpha(0.55)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(2)

        # IQR whiskers
        for arr, pos, color in zip(parts_data, positions, colors):
            q25, q75 = np.percentile(arr, [25, 75])
            iqr = q75 - q25
            lo = max(arr.min(), q25 - 1.5 * iqr)
            hi = min(arr.max(), q75 + 1.5 * iqr)
            ax.vlines(pos, lo, hi, color=color, linewidth=1.5, zorder=3)
            ax.hlines([q25, q75], pos - 0.12, pos + 0.12,
                      color=color, linewidth=2, zorder=4)

    # Significance bracket between the two violins
    y_max = max(
        arr_a.max() if len(arr_a) else 0,
        arr_b.max() if len(arr_b) else 0,
    )
    y_bracket = y_max + 0.04
    ax.annotate(
        "", xy=(2, y_bracket), xytext=(1, y_bracket),
        arrowprops=dict(arrowstyle="-", lw=1.2, color="black"),
    )
    p_txt = f"p={p:.2e}" if not np.isnan(p) else "p=n/a"
    ax.text(
        1.5, y_bracket + 0.01,
        f"{sig_label(p)}  {p_txt}",
        ha="center", va="bottom", fontsize=9,
    )

    _stat_box(ax, arr_a, COLOR_A, x_anchor=0.03)
    _stat_box(ax, arr_b, COLOR_B, x_anchor=0.55)

    handles = [
        Patch(facecolor=COLOR_A, edgecolor="black", label=label_a),
        Patch(facecolor=COLOR_B, edgecolor="black", label=label_b),
    ]
    ax.legend(handles=handles, loc="lower right")

    ax.set_xticks([1, 2])
    ax.set_xticklabels([label_a, label_b])
    ax.set_ylabel("HET / (HET + HOM ALT)")
    ax.set_title("Per-variant Heterozygosity Ratio")
    ax.set_xlim(0.4, 2.6)
    ax.set_ylim(bottom=0)
    ax.margins(y=0.15)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved: {output_path}")

    # Summary stats to TSV
    stats_rows = []
    for label, arr in [(label_a, arr_a), (label_b, arr_b)]:
        stats_rows.append({
            "dataset": label,
            "n": len(arr),
            "mean":   float(np.mean(arr))   if len(arr) else np.nan,
            "median": float(np.median(arr)) if len(arr) else np.nan,
            "Q25":    float(np.percentile(arr, 25)) if len(arr) else np.nan,
            "Q75":    float(np.percentile(arr, 75)) if len(arr) else np.nan,
        })
        print(
            f"  {label}: n={len(arr):,}  median={np.median(arr):.3f}  "
            f"IQR=[{np.percentile(arr, 25):.3f}–{np.percentile(arr, 75):.3f}]"
        )

    stats_df = pd.DataFrame(stats_rows)
    stats_df["U_statistic"] = u_stat
    stats_df["p_value"]     = p
    stats_df["significance"] = sig_label(p)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"Stats saved: {stats_path}")
    print(f"Mann-Whitney U: {u_stat:.1f}  p={p:.4e}  ({sig_label(p)})")


# -- main --

def main():
    parser = argparse.ArgumentParser(
        description="Compare HET / (HET + HOM ALT) distributions between two VCF datasets."
    )
    parser.add_argument("-a", "--vcf-a", required=True, help="First VCF / VCF.gz file")
    parser.add_argument("-b", "--vcf-b", required=True, help="Second VCF / VCF.gz file")
    parser.add_argument("--label-a", default="Dataset A", help="Label for first dataset")
    parser.add_argument("--label-b", default="Dataset B", help="Label for second dataset")
    parser.add_argument(
        "-o", "--out", default="het_ratio_comparison.png",
        help="Output PNG path [default: het_ratio_comparison.png]",
    )
    parser.add_argument(
        "--stats", default="het_ratio_stats.tsv",
        help="Output stats TSV [default: het_ratio_stats.tsv]",
    )
    args = parser.parse_args()

    for path in (args.vcf_a, args.vcf_b):
        if not os.path.exists(path):
            raise FileNotFoundError(f"VCF not found: {path}")

    print(f"Reading {args.vcf_a} ...")
    ratios_a = collect_het_ratios(args.vcf_a)
    print(f"  {len(ratios_a):,} variants with het ratio")

    print(f"Reading {args.vcf_b} ...")
    ratios_b = collect_het_ratios(args.vcf_b)
    print(f"  {len(ratios_b):,} variants with het ratio")

    plot_het_ratio(ratios_a, ratios_b, args.label_a, args.label_b, args.out, args.stats)


if __name__ == "__main__":
    main()
