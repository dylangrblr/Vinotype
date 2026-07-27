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
    python SI_figure_2_variant_quality_metrics_violin.py groupA.vcf.gz groupB.vcf.gz \\
        --label-a Benchmark --label-b Additional \\
        -o comparison.png --stats stats.tsv
"""

from __future__ import annotations

import argparse
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import mannwhitneyu

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def open_vcf(path: str):
    """Open plain or gzipped VCF."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str: str) -> dict:
    info: dict = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v
        else:
            info[field] = True
    return info


def _float(val: str):
    """Return float or None for missing/empty values."""
    v = str(val).split(",")[0].strip()
    if v in (".", "", "NA"):
        return None
    try:
        return float(v)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Per-variant extraction
# ---------------------------------------------------------------------------

def extract_record(parts: list) -> dict:
    """
    Return a dict with per-site keys (AC, AN, AF, MAF) and per-sample
    aggregated keys (DP_wt, DP_alt, GQ_wt, GQ_alt, VAF_alt).
    _wt = mean across 0/0 samples; _alt = mean across alt-carrying samples.

    FORMAT fields are accessed by pre-computed integer indices
    """
    if len(parts) < 10:   # need FORMAT + at least one sample
        return {}

    info = parse_info(parts[7])
    result: dict = {}

    # --- Per-site metrics --------------------------------------------------
    result["AC"] = _float(info["AC"]) if "AC" in info else None
    result["AN"] = _float(info["AN"]) if "AN" in info else None

    af = _float(info["AF"]) if "AF" in info else None
    if af is None and result["AC"] is not None and result["AN"]:
        af = result["AC"] / result["AN"]
    result["AF"] = af

    maf = _float(info["MAF"]) if "MAF" in info else None
    if maf is None and af is not None:
        maf = min(af, 1.0 - af)
    result["MAF"] = maf

    # --- Pre-compute FORMAT field indices once for all samples -------------
    fmt_keys = parts[8].split(":")
    fki      = {k: i for i, k in enumerate(fmt_keys)}
    gt_i     = fki.get("GT",  -1)
    dp_i     = fki.get("DP",  -1)
    gq_i     = fki.get("GQ",  -1)
    ad_i     = fki.get("AD",  -1)
    vaf_i    = fki.get("VAF", -1)
    ao_i     = fki.get("AO",  -1)
    ro_i     = fki.get("RO",  -1)

    if gt_i < 0:
        return result   # no GT field — cannot stratify

    # --- Aggregate per-sample metrics across all samples -------------------
    wt_dp,  alt_dp  = [], []
    wt_gq,  alt_gq  = [], []
    alt_vaf         = []

    for sample_str in parts[9:]:
        fields = sample_str.split(":")
        nf = len(fields)

        if gt_i >= nf:
            continue
        gt_raw = fields[gt_i]
        if gt_raw in (".", "./.", ".|."):
            continue
        gt    = gt_raw.replace("|", "/")
        valid = [a for a in gt.split("/") if a != "."]
        if not valid:
            continue
        is_alt = any(a != "0" for a in valid)

        # DP
        if 0 <= dp_i < nf:
            dp = _float(fields[dp_i])
            if dp is not None:
                (alt_dp if is_alt else wt_dp).append(dp)

        # GQ
        if 0 <= gq_i < nf:
            gq = _float(fields[gq_i])
            if gq is not None:
                (alt_gq if is_alt else wt_gq).append(gq)

        # VAF —> only alt carrying samples
        if is_alt:
            vaf = None
            # 1. Explicit VAF tag
            if 0 <= vaf_i < nf:
                vaf = _float(fields[vaf_i])
            # 2. AD-derived
            if vaf is None and 0 <= ad_i < nf:
                raw = fields[ad_i]
                if raw not in (".", ""):
                    try:
                        ad = [int(x) for x in raw.split(",")]
                        total = sum(ad)
                        if total > 0 and len(ad) >= 2:
                            vaf = sum(ad[1:]) / total
                    except ValueError:
                        pass
            # 3. AO/RO-derived
            if vaf is None and 0 <= ao_i < nf and 0 <= ro_i < nf:
                try:
                    ao = sum(int(x) for x in fields[ao_i].split(","))
                    ro = int(fields[ro_i])
                    total = ao + ro
                    if total > 0:
                        vaf = ao / total
                except ValueError:
                    pass
            # 4. GT dosage fallback
            if vaf is None:
                vaf = sum(1 for a in valid if a != "0") / len(valid)
            alt_vaf.append(vaf)

    # Store per-locus means (None when no samples contributed)
    result["DP_wt"]   = float(np.mean(wt_dp))   if wt_dp   else None
    result["DP_alt"]  = float(np.mean(alt_dp))  if alt_dp  else None
    result["GQ_wt"]   = float(np.mean(wt_gq))   if wt_gq   else None
    result["GQ_alt"]  = float(np.mean(alt_gq))  if alt_gq  else None
    result["VAF_alt"] = float(np.mean(alt_vaf)) if alt_vaf else None

    return result


# ---------------------------------------------------------------------------
# VCF reader
# ---------------------------------------------------------------------------

# Ordered list
METRICS = [
    "AC", "AN", "AF", "MAF",
    "DP_wt", "DP_alt",
    "GQ_wt", "GQ_alt",
    "VAF_alt",
]

def read_vcf(path: str) -> dict[str, list]:
    """Return {metric: [values]} for a single VCF file."""
    data: dict[str, list] = {m: [] for m in METRICS}

    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            rec = extract_record(parts)
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
    "AC":      ("Allele Count",                           "royalblue"),
    "AN":      ("Allele Number",                          "steelblue"),
    "AF":      ("Allele Frequency",                       "seagreen"),
    "MAF":     ("Minor Allele Frequency",                 "mediumseagreen"),
    "DP_wt":   ("Read Depth — WT samples (0/0)",          "#e07b54"),
    "DP_alt":  ("Read Depth — alt-carrying samples",      "tomato"),
    "GQ_wt":   ("Genotype Quality — WT samples (0/0)",    "#9b82c4"),
    "GQ_alt":  ("Genotype Quality — alt-carrying samples", "mediumpurple"),
    "VAF_alt": ("Variant Allele Frequency — alt-carrying samples", "darkorange"),
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
        fontsize=7,
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
        if len(arr_a):
            parts_data.append(arr_a)
            parts_pos.append(1)
            parts_colors.append(COLOR_A)
            parts_labels.append(label_a)
        if len(arr_b):
            parts_data.append(arr_b)
            parts_pos.append(2)
            parts_colors.append(COLOR_B)
            parts_labels.append(label_b)

        if parts_data:
            vp = ax.violinplot(
                parts_data,
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
                fontsize=9, fontweight="bold",
            )
        else:
            ax.set_title(desc, fontsize=9, fontweight="bold")

        ax.set_ylabel(metric, fontsize=9)
        ax.set_xticks(parts_pos)
        ax.set_xticklabels(parts_labels, fontsize=9)

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
                fontsize=7,
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
