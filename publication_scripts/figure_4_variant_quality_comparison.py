#!/usr/bin/env python3
"""
Compare per-variant quality metrics (DP, GQ, MQ, MAF) between two VCF datasets
and visualise distributions as paired violin plots with Mann-Whitney U tests.

  DP  — read depth, averaged across heterozygous samples per locus (FORMAT/DP)
  GQ  — genotype quality, averaged across heterozygous samples per locus (FORMAT/GQ)
  MQ  — RMS mapping quality (INFO/MQ)
  MAF — minor allele frequency (INFO/MAF, or derived from INFO/AF)

Usage:
    python3 figure_4_variant_quality_comparison.py \\
        -a datasetA.vcf.gz -b datasetB.vcf.gz \\
        --label-a DatasetA --label-b DatasetB \\
        -o quality_comparison.png --stats quality_stats.tsv
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

# -- metric extraction (adapted from figure_6_plot_variant_allele_freq.py) --

# cyvcf2 sentinel for missing integer FORMAT fields
_MISSING_INT = -2147483648


def collect_metrics(vcf_path):
    """Return {metric: [values]} by iterating over every variant in a VCF."""
    data = {"DP": [], "GQ": [], "MQ": [], "MAF": []}

    for variant in VCF(vcf_path):

        # MQ — INFO field
        mq = variant.INFO.get("MQ")
        if mq is not None:
            try:
                data["MQ"].append(float(mq))
            except (TypeError, ValueError):
                pass

        # MAF — prefer INFO/MAF, fall back to min(AF, 1−AF)
        maf = variant.INFO.get("MAF")
        if maf is None:
            af_raw = variant.INFO.get("AF")
            if af_raw is not None:
                try:
                    af = float(str(af_raw).split(",")[0])
                    maf = min(af, 1.0 - af)
                except (TypeError, ValueError):
                    maf = None
        if maf is not None:
            try:
                data["MAF"].append(float(maf))
            except (TypeError, ValueError):
                pass

        # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
        het_mask = variant.gt_types == 1

        # DP — FORMAT field, mean across heterozygous samples only
        dp_fmt = variant.format("DP")
        if dp_fmt is not None:
            vals = [
                int(dp_fmt[i][0])
                for i in range(len(dp_fmt))
                if het_mask[i] and dp_fmt[i][0] != _MISSING_INT and dp_fmt[i][0] > 0
            ]
            if vals:
                data["DP"].append(float(np.mean(vals)))

        # GQ — FORMAT field, mean across heterozygous samples only
        gq_fmt = variant.format("GQ")
        if gq_fmt is not None:
            vals = [
                int(gq_fmt[i][0])
                for i in range(len(gq_fmt))
                if het_mask[i] and gq_fmt[i][0] != _MISSING_INT and gq_fmt[i][0] >= 0
            ]
            if vals:
                data["GQ"].append(float(np.mean(vals)))

    return data


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

METRIC_META = {
    "DP":  "Mean Read Depth (per locus)",
    "GQ":  "Mean Genotype Quality (per locus)",
    "MQ":  "RMS Mapping Quality",
    "MAF": "Minor Allele Frequency",
}


def _stat_box(ax, arr, color, x_anchor):
    """Summary stats inset (borrowed from SI_figure_2_variant_quality_metrics_violin.py)."""
    if len(arr) == 0:
        return
    txt = (
        f"n={len(arr):,}\n"
        f"med={np.median(arr):.3g}\n"
        f"IQR [{np.percentile(arr, 25):.3g}–{np.percentile(arr, 75):.3g}]"
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


def plot_quality(data_a, data_b, label_a, label_b, output_path, stats_path):
    metrics   = [m for m in ["DP", "GQ", "MQ", "MAF"] if data_a[m] or data_b[m]]
    n_metrics = len(metrics)

    if n_metrics == 0:
        print("No plottable metrics found in either VCF.")
        return

    ncols = 2
    nrows = (n_metrics + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 5 * nrows), squeeze=False)

    stats_rows = []

    for idx, metric in enumerate(metrics):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        arr_a = np.array(data_a[metric]) if data_a[metric] else np.array([])
        arr_b = np.array(data_b[metric]) if data_b[metric] else np.array([])

        parts_data, positions, colors = [], [], []
        if len(arr_a):
            parts_data.append(arr_a)
            positions.append(1)
            colors.append(COLOR_A)
        if len(arr_b):
            parts_data.append(arr_b)
            positions.append(2)
            colors.append(COLOR_B)

        # Clip DP to 1st–99th percentile so outliers don't distort KDE shape
        if metric == "DP" and parts_data:
            _comb = np.concatenate(parts_data)
            lo, hi = np.percentile(_comb, [1, 99])
            parts_data = [
                arr[(arr >= lo) & (arr <= hi)] if len(arr) else arr
                for arr in parts_data
            ]

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

            # IQR whiskers (full arrays, not clipped)
            for arr_raw, pos, color in zip(
                [arr_a, arr_b] if len(arr_a) and len(arr_b) else parts_data,
                positions, colors,
            ):
                if len(arr_raw) == 0:
                    continue
                q25, q75 = np.percentile(arr_raw, [25, 75])
                iqr = q75 - q25
                lo_w = max(arr_raw.min(), q25 - 1.5 * iqr)
                hi_w = min(arr_raw.max(), q75 + 1.5 * iqr)
                ax.vlines(pos, lo_w, hi_w, color=color, linewidth=1.5, zorder=3)
                ax.hlines([q25, q75], pos - 0.12, pos + 0.12,
                          color=color, linewidth=2, zorder=4)

        # Significance bracket
        u_stat, p = mw_test(arr_a.tolist(), arr_b.tolist())
        if len(arr_a) and len(arr_b):
            y_top = np.percentile(np.concatenate([arr_a, arr_b]), 99) * 1.06
            ax.annotate(
                "", xy=(2, y_top), xytext=(1, y_top),
                arrowprops=dict(arrowstyle="-", lw=1.2, color="black"),
            )
            p_txt = f"p={p:.2e}" if not np.isnan(p) else "p=n/a"
            ax.text(1.5, y_top * 1.01, f"{sig_label(p)}  {p_txt}",
                    ha="center", va="bottom", fontsize=8)

        _stat_box(ax, arr_a, COLOR_A, x_anchor=0.03)
        _stat_box(ax, arr_b, COLOR_B, x_anchor=0.55)

        if not len(arr_a):
            ax.text(1, 0.5, "no data", transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=9,
                    color=COLOR_A, style="italic", alpha=0.7)
        if not len(arr_b):
            ax.text(2, 0.5, "no data", transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=9,
                    color=COLOR_B, style="italic", alpha=0.7)

        ax.set_xticks([1, 2])
        ax.set_xticklabels([label_a, label_b])
        ax.set_ylabel(METRIC_META[metric])
        ax.set_title(METRIC_META[metric], fontsize=9, fontweight="bold")
        ax.set_xlim(0.4, 2.6)

        stats_rows.append({
            "metric":       metric,
            f"n_{label_a}":      len(arr_a),
            f"median_{label_a}": float(np.median(arr_a)) if len(arr_a) else np.nan,
            f"n_{label_b}":      len(arr_b),
            f"median_{label_b}": float(np.median(arr_b)) if len(arr_b) else np.nan,
            "U_statistic":  u_stat,
            "p_value":      p,
            "significance": sig_label(p),
        })

    # Hide unused subplot slots
    for idx in range(n_metrics, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    handles = [
        Patch(facecolor=COLOR_A, edgecolor="black", label=label_a),
        Patch(facecolor=COLOR_B, edgecolor="black", label=label_b),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved: {output_path}")

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"Stats saved: {stats_path}")


# -- main --

def main():
    parser = argparse.ArgumentParser(
        description="Compare DP, GQ, MQ, and MAF distributions between two VCF datasets."
    )
    parser.add_argument("-a", "--vcf-a", required=True, help="First VCF / VCF.gz file")
    parser.add_argument("-b", "--vcf-b", required=True, help="Second VCF / VCF.gz file")
    parser.add_argument("--label-a", default="Dataset A", help="Label for first dataset")
    parser.add_argument("--label-b", default="Dataset B", help="Label for second dataset")
    parser.add_argument(
        "-o", "--out", default="quality_comparison.png",
        help="Output PNG path [default: quality_comparison.png]",
    )
    parser.add_argument(
        "--stats", default="quality_stats.tsv",
        help="Output stats TSV [default: quality_stats.tsv]",
    )
    args = parser.parse_args()

    for path in (args.vcf_a, args.vcf_b):
        if not os.path.exists(path):
            raise FileNotFoundError(f"VCF not found: {path}")

    print(f"Reading {args.vcf_a} ...")
    data_a = collect_metrics(args.vcf_a)
    print(f"  DP={len(data_a['DP']):,}  GQ={len(data_a['GQ']):,}  "
          f"MQ={len(data_a['MQ']):,}  MAF={len(data_a['MAF']):,} variants")

    print(f"Reading {args.vcf_b} ...")
    data_b = collect_metrics(args.vcf_b)
    print(f"  DP={len(data_b['DP']):,}  GQ={len(data_b['GQ']):,}  "
          f"MQ={len(data_b['MQ']):,}  MAF={len(data_b['MAF']):,} variants")

    plot_quality(data_a, data_b, args.label_a, args.label_b, args.out, args.stats)


if __name__ == "__main__":
    main()
