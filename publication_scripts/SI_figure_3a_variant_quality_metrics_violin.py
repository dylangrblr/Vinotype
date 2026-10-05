#!/usr/bin/env python3
"""
CCompare variant calling metrics between two datasets. Group A is treated as the 
benchmarking dataset (variants identified by both benchmarking and Vinotype), and 
Group B as the additional (variants identified by Vinotype only). The two groups 
are assumed to be from the same dataset. For each metric the two distributions are 
compared with a two-tailed Mann-Whitney U test and visualised with paired violin plots.

Per-site metrics  : AC, AN, AF, MAF
Per-sample metrics: DP_wt, DP_alt, GQ_wt, GQ_alt, VAF_alt
                    Averaged across ALL samples at each locus, split by
                    genotype class:
                      _wt  = wild-type (0/0) samples only
                      _alt = alt-carrying (0/1, 1/1, …) samples only

Usage:
    python SI_figure_2_variant_quality_metrics_violin.py \\
        groupA.vcf.gz \\
        groupB.vcf.gz \\
        --label-a Benchmark \\
        --label-b Additional \\
        -o comparison.png \\
        --stats stats.tsv
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from cyvcf2 import VCF
from scipy.stats import mannwhitneyu

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _first_scalar(val):
    """Take the first element of an INFO value (handles arrays/scalars)."""
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, (tuple, list, np.ndarray)):
        return float(val[0]) if len(val) else None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Per-variant extraction (cyvcf2, vectorized)
# ---------------------------------------------------------------------------

def extract_record(var) -> dict:
    """
    Return a dict with per-site keys (AC, AN, AF, MAF) and per-sample
    aggregated keys split by genotype class (_wt / _het / _hom).
    Uses cyvcf2 NumPy accessors — no per-sample Python string parsing.
    """
    info = var.INFO
    result: dict = {}

    # --- Per-site metrics --------------------------------------------------
    result["AC"] = _first_scalar(info.get("AC"))
    result["AN"] = _first_scalar(info.get("AN"))

    af = _first_scalar(info.get("AF"))
    if af is None and result["AC"] is not None and result["AN"]:
        af = result["AC"] / result["AN"]
    result["AF"] = af

    maf = _first_scalar(info.get("MAF"))
    if maf is None and af is not None:
        maf = min(af, 1.0 - af)
    result["MAF"] = maf

    result["MQ"]  = _float(info["MQ"])  if "MQ"  in info else None
    result["SOR"] = _float(info["SOR"]) if "SOR" in info else None

    result["MQ"]  = _first_scalar(info.get("MQ"))
    result["SOR"] = _first_scalar(info.get("SOR"))

    # --- Per-sample metrics (vectorized) -----------------------------------
    gt_types = var.gt_types           # 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
    known = gt_types != 2
    if not known.any():
        return result                 # no GT field or all missing

    wt  = known & (gt_types == 0)
    het = known & (gt_types == 1)
    hom = known & (gt_types == 3)

    het_count = int(het.sum())
    hom_count = int(hom.sum())
    _total_alt = het_count + hom_count
    result["HetRatio"] = het_count / _total_alt if _total_alt > 0 else None

    # DP / GQ: missing values come back as -1 (DP int) / -1.0 (GQ float).
    # Guard against undeclared tags (KeyError from cyvcf2 header lookup).
    n_samples = gt_types.shape[0]
    try:
        dp = np.asarray(var.gt_depths, dtype=float)
    except KeyError:
        dp = np.full(n_samples, np.nan)
    dp = np.where(dp >= 0, dp, np.nan)
    try:
        gq = np.asarray(var.gt_quals, dtype=float)
    except KeyError:
        gq = np.full(n_samples, np.nan)
    gq = np.where(gq >= 0, gq, np.nan)

    def _masked_mean(arr: np.ndarray, mask: np.ndarray):
        sel = arr[mask]
        sel = sel[~np.isnan(sel)]
        return float(sel.mean()) if sel.size else None

    result["DP_wt"]  = _masked_mean(dp, wt)
    result["DP_het"] = _masked_mean(dp, het)
    result["DP_hom"] = _masked_mean(dp, hom)
    result["GQ_wt"]  = _masked_mean(gq, wt)
    result["GQ_het"] = _masked_mean(gq, het)
    result["GQ_hom"] = _masked_mean(gq, hom)

    # VAF — alt-carrying samples only, split HET / HOM ALT.
    # Source: explicit VAF tag > AD-derived > AO/RO-derived. Samples with none
    # of these stay NaN and are excluded from the means.
    vaf_arr = np.full(gt_types.shape[0], np.nan)

    def _fmt(name):
        """var.format() but returns None when the tag is absent/undeclared."""
        try:
            return var.format(name)
        except KeyError:
            return None

    tag = _fmt("VAF")                 # explicit VAF FORMAT tag
    if tag is not None:
        vaf_arr = np.asarray(tag, dtype=float).reshape(-1)

    if np.isnan(vaf_arr).all():       # fall back to AD
        ad = _fmt("AD")               # (n_samples, n_alleles)
        if ad is not None:
            ad = np.asarray(ad, dtype=float)
            ad = np.where(ad >= 0, ad, np.nan)
            if ad.ndim > 1 and ad.shape[1] >= 2:
                total = np.nansum(ad, axis=1)
                alt   = np.nansum(ad[:, 1:], axis=1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    vaf_arr = np.where(total > 0, alt / total, np.nan)

    if np.isnan(vaf_arr).all():       # fall back to AO/RO
        ao = _fmt("AO")
        ro = _fmt("RO")
        if ao is not None and ro is not None:
            ao = np.asarray(ao, dtype=float)
            ao = np.where(ao >= 0, ao, np.nan)
            ao_sum = np.nansum(ao, axis=1) if ao.ndim > 1 else ao
            ro = np.asarray(ro, dtype=float).reshape(-1)
            ro = np.where(ro >= 0, ro, np.nan)
            total = ao_sum + ro
            with np.errstate(invalid="ignore", divide="ignore"):
                vaf_arr = np.where(total > 0, ao_sum / total, np.nan)

    # Samples with a VAF value still missing (no VAF/AD/AO/RO) are excluded
    # from the per-locus means below via NaN filtering.

    het_vaf = vaf_arr[het]
    het_vaf = het_vaf[~np.isnan(het_vaf)]
    hom_vaf = vaf_arr[hom]
    hom_vaf = hom_vaf[~np.isnan(hom_vaf)]

    result["VAF_het"] = float(het_vaf.mean()) if het_vaf.size else None
    result["VAF_hom"] = float(hom_vaf.mean()) if hom_vaf.size else None

    return result


# ---------------------------------------------------------------------------
# VCF reader
# ---------------------------------------------------------------------------

# Ordered list
METRICS = [
    # "AC",
    # "AN",
    # "AF",
    "MAF",
    "MQ",
    "SOR",
    "HetRatio",
    "DP_wt", "DP_het", "DP_hom",
    "GQ_wt", "GQ_het", "GQ_hom",
    "VAF_het", "VAF_hom",
]

def read_vcf(path: str) -> dict[str, list]:
    """Return {metric: [values]} for a single VCF file (plain or bgzipped)."""
    data: dict[str, list] = {m: [] for m in METRICS}

    for var in VCF(path):             # cyvcf2 handles plain/gzip transparently
        rec = extract_record(var)
        for m in METRICS:
            val = rec.get(m)
            if val is not None:
                data[m].append(val)

    return data


# ---------------------------------------------------------------------------
# Mann-Whitney U (two-tailed)
# ---------------------------------------------------------------------------

def mw_test(a: list, b: list) -> tuple[float, float]:
    """Return (U-statistic, two-tailed p-value). NaN when insufficient data."""
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = mannwhitneyu(a, b, alternative="two-sided", use_continuity=True)
    return res.statistic, res.pvalue


# ---------------------------------------------------------------------------
# Significance labels
# ---------------------------------------------------------------------------

def sig_label(p: float) -> str:
    if np.isnan(p):
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

METRIC_META = {
    "AC":       ("Allele Count",                                        "royalblue"),
    "AN":       ("Allele Number",                                       "steelblue"),
    "AF":       ("Allele Frequency",                                    "seagreen"),
    "MAF":      ("Minor Allele Frequency",                              "mediumseagreen"),
    "MQ":       ("RMS Mapping Quality",                                 "darkslategray"),
    "SOR":      ("Strand Odds Ratio",                                   "firebrick"),
    "HetRatio": ("Heterozygosity Ratio [HET / (HET + HOM ALT)]",       "goldenrod"),
    "DP_wt":    ("Read Depth — WT samples (0/0)",                       "#e07b54"),
    "DP_het":   ("Read Depth — heterozygous samples",                   "tomato"),
    "DP_hom":   ("Read Depth — homozygous alt samples",                 "darkred"),
    "GQ_wt":    ("Genotype Quality — WT samples (0/0)",                 "#9b82c4"),
    "GQ_het":   ("Genotype Quality — heterozygous samples",             "mediumpurple"),
    "GQ_hom":   ("Genotype Quality — homozygous alt samples",           "purple"),
    "VAF_het":  ("Variant Allele Frequency — heterozygous samples",     "darkorange"),
    "VAF_hom":  ("Variant Allele Frequency — homozygous alt samples",   "orangered"),
}

# Colours for the two groups
COLOR_A = "#4C72B0"   # blue  → benchmark
COLOR_B = "#DD8452"   # orange → additional


def _stat_box(ax, arr, color):
    """Small summary text box."""
    if len(arr) == 0:
        return
    txt = (
        f"n={len(arr):,}\n"
        f"med={np.median(arr):.3g}\n"
        f"IQR [{np.percentile(arr,25):.3g}–{np.percentile(arr,75):.3g}]"
    )
    ax.text(
        0.03, 0.97, txt,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        color=color,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor=color, alpha=0.7),
    )


def plot_comparison(
    data_a: dict,
    data_b: dict,
    label_a: str,
    label_b: str,
    stats_df: pd.DataFrame,
    output_path: str,
):
    available = [m for m in METRICS if data_a[m] or data_b[m]]
    n = len(available)

    if n == 0:
        print("No plottable metrics found.")
        return

    ncols = 2
    nrows = (n + 1) // 2
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(14, 4.5 * nrows),
        squeeze=False,
    )

    for idx, metric in enumerate(available):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        desc, _ = METRIC_META[metric]

        arr_a = np.array(data_a[metric]) if data_a[metric] else np.array([])
        arr_b = np.array(data_b[metric]) if data_b[metric] else np.array([])

        # Build violin data list (skip empty arrays)
        parts_data, parts_pos, parts_colors, parts_labels = [], [], [], []
        violin_positions = (1.15, 1.85)
        if len(arr_a):
            parts_data.append(arr_a)
            parts_pos.append(violin_positions[0])
            parts_colors.append(COLOR_A)
            parts_labels.append(label_a)
        if len(arr_b):
            parts_data.append(arr_b)
            parts_pos.append(violin_positions[1])
            parts_colors.append(COLOR_B)
            parts_labels.append(label_b)

        # For DP metrics, filter to the 1st–99th percentile before KDE so
        # the violin shape is not distorted by outliers outside the visible range.
        if metric.startswith("DP_") and parts_data:
            _comb = np.concatenate(parts_data)
            _lo_dp, _hi_dp = np.percentile(_comb, [1, 99])
            violin_data = []
            for arr in parts_data:
                clipped = arr[(arr >= _lo_dp) & (arr <= _hi_dp)]
                violin_data.append(clipped if len(clipped) > 0 else arr)
        else:
            violin_data = parts_data

        if parts_data:
            vp = ax.violinplot(
                violin_data,
                positions=parts_pos,
                widths=0.6,
                showmedians=True,
                showextrema=False,
            )
            for body, color in zip(vp["bodies"], parts_colors):
                body.set_facecolor(color)
                body.set_alpha(0.55)
            vp["cmedians"].set_color("black")
            vp["cmedians"].set_linewidth(2)

            # Overlay box-plot whiskers
            for arr, pos, color in zip(parts_data, parts_pos, parts_colors):
                q25, q75 = np.percentile(arr, [25, 75])
                iqr = q75 - q25
                lo = max(arr.min(), q25 - 1.5 * iqr)
                hi = min(arr.max(), q75 + 1.5 * iqr)
                ax.vlines(pos, lo, hi, color=color, linewidth=1.5, zorder=3)
                ax.hlines([q25, q75], pos - 0.12, pos + 0.12,
                          color=color, linewidth=2, zorder=4)

        # Significance annotation
        row_s = stats_df[stats_df["Metric"] == metric]
        if not row_s.empty:
            p = float(row_s["p_value"].values[0])
            u = float(row_s["U_statistic"].values[0])
            label_sig = sig_label(p)
            p_txt = f"p={p:.2e}" if not np.isnan(p) else "p=n/a"
            ax.set_title(
                f"{desc}  [{label_sig}  {p_txt}]",
                fontsize=11, fontweight="bold",
            )
        else:
            ax.set_title(desc, fontsize=11, fontweight="bold")

        ax.set_ylabel(metric, fontsize=11)
        ax.set_xticks(violin_positions)
        ax.set_xticklabels([label_a, label_b], fontsize=11)
        ax.tick_params(axis="y", labelsize=10)
        # Mark positions where a dataset has no data for this metric
        for _nd_arr, _nd_pos, _nd_col in [
            (arr_a, violin_positions[0], COLOR_A),
            (arr_b, violin_positions[1], COLOR_B),
        ]:
            if not len(_nd_arr):
                ax.text(
                    _nd_pos, 0.5, "no data",
                    transform=ax.get_xaxis_transform(),
                    ha="center", va="center",
                    fontsize=9, color=_nd_col, style="italic", alpha=0.7,
                )

        # Stats boxes
        if len(arr_a):
            _stat_box(ax, arr_a, COLOR_A)
        if len(arr_b) and len(arr_a):
            # Second box on right side
            txt_b = (
                f"n={len(arr_b):,}\n"
                f"med={np.median(arr_b):.3g}\n"
                f"IQR [{np.percentile(arr_b,25):.3g}–{np.percentile(arr_b,75):.3g}]"
            )
            ax.text(
                0.97, 0.97, txt_b,
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment="top",
                horizontalalignment="right",
                color=COLOR_B,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=COLOR_B, alpha=0.7),
            )
        elif len(arr_b):
            _stat_box(ax, arr_b, COLOR_B)

        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f"{x:,.3g}")
        )

        # Clip y-axis for DP metrics
        if metric.startswith("DP_") and parts_data:
            combined = np.concatenate(parts_data)
            lo = np.percentile(combined, 1)
            hi = np.percentile(combined, 99)
            pad = max((hi - lo) * 0.20, 1.0)
            ax.set_ylim(max(0.0, lo - pad), hi + pad)

    # Hide unused panels
    for idx in range(len(available), nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle(
        f"Variant metric comparison\n"
        f"{label_a} (benchmark) vs {label_b} (additional)  —  "
        f"two-tailed Mann-Whitney U\n"
        f"Sample-level metrics (DP, GQ, VAF) are per-locus means "
        f"across all samples, split by genotype class (_wt = 0/0; _alt = alt-carrying)",
        fontsize=11, fontweight="bold", y=1.02,
    )
    # Legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=COLOR_A, alpha=0.7, label=label_a),
        Patch(facecolor=COLOR_B, alpha=0.7, label=label_b),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=10,
        framealpha=0.9,
        bbox_to_anchor=(1.0, 1.0),
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved  → {output_path}")


# ---------------------------------------------------------------------------
# Console summary table
# ---------------------------------------------------------------------------

def print_stats_table(df: pd.DataFrame):
    col_w = {"Metric": 6, "n_A": 9, "n_B": 9,
             "median_A": 11, "median_B": 11,
             "U_statistic": 13, "p_value": 12, "sig": 4}
    header = "  ".join(c.ljust(col_w[c]) for c in df.columns)
    print("\n" + "=" * len(header))
    print("Mann-Whitney U results  (two-tailed)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for _, r in df.iterrows():
        def fmt(c):
            v = r[c]
            if c in ("n_A", "n_B"):
                return str(int(v)).ljust(col_w[c])
            if c in ("median_A", "median_B", "U_statistic"):
                return (f"{v:.4g}" if not np.isnan(v) else "n/a").ljust(col_w[c])
            if c == "p_value":
                return (f"{v:.3e}" if not np.isnan(v) else "n/a").ljust(col_w[c])
            return str(v).ljust(col_w[c])
        print("  ".join(fmt(c) for c in df.columns))
    print("=" * len(header))
    print("Significance: *** p<0.001  ** p<0.01  * p<0.05  ns p≥0.05\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare variant quality metrics (AC, AN, AF, MAF, AD, DP, GQ, VAF) "
            "between two VCF files using the two-tailed Mann-Whitney U test."
        )
    )
    parser.add_argument("vcf_a", help="Group A VCF (benchmark dataset)")
    parser.add_argument("vcf_b", help="Group B VCF (additional dataset)")
    parser.add_argument(
        "--label-a", default="Group A (Benchmark)",
        help="Display label for group A (default: 'Group A (Benchmark)')",
    )
    parser.add_argument(
        "--label-b", default="Group B (Additional)",
        help="Display label for group B (default: 'Group B (Additional)')",
    )
    parser.add_argument(
        "-o", "--output", default="vcf_comparison.png",
        help="Output plot file (default: vcf_comparison.png)",
    )
    parser.add_argument(
        "--no-numpy-warn", action="store_true",
        help="Suppress numpy mean-of-empty-slice warnings.",
    )
    parser.add_argument(
        "--stats", default=None,
        help="Optional TSV file to write the statistics table to.",
    )
    args = parser.parse_args()

    if args.no_numpy_warn:
        warnings.filterwarnings("ignore", category=RuntimeWarning)

    print(f"Reading group A: {args.vcf_a}")
    data_a = read_vcf(args.vcf_a)
    print(f"  Variants parsed: {max(len(v) for v in data_a.values()):,}")

    print(f"\nReading group B: {args.vcf_b}")
    data_b = read_vcf(args.vcf_b)
    print(f"  Variants parsed: {max(len(v) for v in data_b.values()):,}")

    # -----------------------------------------------------------------------
    # Build statistics table
    # -----------------------------------------------------------------------
    rows = []
    for metric in METRICS:
        a = data_a[metric]
        b = data_b[metric]
        u_stat, p_val = mw_test(a, b)
        rows.append({
            "Metric":      metric,
            "n_A":         len(a),
            "n_B":         len(b),
            "median_A":    np.median(a) if a else float("nan"),
            "median_B":    np.median(b) if b else float("nan"),
            "U_statistic": u_stat,
            "p_value":     p_val,
            "sig":         sig_label(p_val),
        })

    stats_df = pd.DataFrame(rows)
    print_stats_table(stats_df)

    if args.stats:
        stats_df.to_csv(args.stats, sep="\t", index=False)
        print(f"Stats table saved → {args.stats}")

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    plot_comparison(
        data_a, data_b,
        args.label_a, args.label_b,
        stats_df,
        args.output,
    )


if __name__ == "__main__":
    main()
