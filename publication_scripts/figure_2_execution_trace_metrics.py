#!/usr/bin/env python3
"""
Summarise Nextflow resource usage from a pipeline execution trace and plot
key metrics (wall-clock time, CPU utilisation, peak memory, I/O volume) as
box plots, one box per individual process, coloured/grouped by broad
pipeline stage. Values are plotted on their real (untransformed) scale.

Usage:
    python3 figure_2_execution_trace_metrics.py \\
        -i execution_trace.txt \\
        -o execution_trace_metrics.png \\
        --stats execution_trace_stats.tsv
"""

import argparse
import colorsys
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# -- broad stage grouping ------------------------------------------------

# Order controls left-to-right placement on the plots (pipeline order).
STAGE_ORDER = [
    "QUALITY_CONTROL",
    "REFERENCE_ALIGNMENT",
    "VARIANT_CALLING",
    "JOINT_GENOTYPING",
    "VARIANT_FILTERING",
    "VARIANT_EFFECT_PREDICTION",
    "MULTIQC",
]

STAGE_LABELS = {
    "QUALITY_CONTROL": "Quality control",
    "REFERENCE_ALIGNMENT": "Reference alignment",
    "VARIANT_CALLING": "Variant calling",
    "JOINT_GENOTYPING": "Joint genotyping",
    "VARIANT_FILTERING": "Variant filtering",
    "VARIANT_EFFECT_PREDICTION": "Variant effect prediction",
    "MULTIQC": "MultiQC",
}

# Per-stage process call order, taken from each subworkflow's main.nf (not alphabetical/
# guessed), so panel x-axes reflect the true pipeline execution sequence.
PROCESS_ORDER = {
    "QUALITY_CONTROL": ["FASTQC_PRE", "FASTP"],
    "REFERENCE_ALIGNMENT": [
        "SAMTOOLS_FAIDX", "SAMTOOLS_DICT", "GENOZIP_REFERENCE",
        "BWAMEM2_MEM", "PICARD_SORTSAM", "PICARD_ADDORREPLACEREADGROUPS",
        "PICARD_MARKDUPLICATES", "SAMTOOLS_FLAGSTAT_PRE", "SAMTOOLS_VIEW",
        "SAMTOOLS_FLAGSTAT_POST", "SAMTOOLS_STATS", "GENOZIP_COMPRESS__BAM",
    ],
    "VARIANT_CALLING": ["GATK4_HAPLOTYPECALLER", "GENOZIP_COMPRESS__GVCF", "BCFTOOLS_SPLITGVCFBYCHROMOSOME"],
    "JOINT_GENOTYPING": ["GATK4_REBLOCKGVCF", "GATK4_GENOMICSDBIMPORT", "GATK4_GENOTYPEGVCFS"],
    "VARIANT_FILTERING": [
        "GATK4_SELECTVARIANTS__SNP", "GATK4_SELECTVARIANTS__INDEL", "GATK4_VARIANTSTOTABLE",
        "BCFTOOLS_INDEX", "BCFTOOLS_QUERY", "UTILS_FINALFILTER",
    ],
    "VARIANT_EFFECT_PREDICTION": ["GATK4_MERGEVCFS", "SNPEFF_SNPEFF", "SNPSIFT_FILTER"],
    "MULTIQC": ["MULTIQC"],
}


def broad_stage(name):
    """Map a fully-qualified Nextflow process name to a broad pipeline stage."""
    # e.g. "DYLANGRBLR:VINOTYPE:REFERENCE_ALIGNMENT:INDEXES:SAMTOOLS_FAIDX (T2T_ref.fasta)"
    parts = name.split(" (")[0].split(":")
    for part in parts:
        if part in STAGE_ORDER:
            return part
    return None


def extract_process(name):
    """Return the specific Nextflow process (module) name, e.g. 'BWAMEM2_MEM'."""
    return name.split(" (")[0].split(":")[-1]


# -- unit parsing ----------------------------------------------------------

_DUR_RE = re.compile(r"(?:(\d+)h)?\s*(?:(\d+)m)?\s*(?:([\d.]+)s)?")

_SIZE_UNITS = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


def parse_duration(value):
    """Convert Nextflow duration strings (e.g. '6m 53s', '1h 2m 3s') to seconds."""
    if pd.isna(value) or value in ("-", ""):
        return np.nan
    match = _DUR_RE.fullmatch(value.strip())
    if not match or not any(match.groups()):
        return np.nan
    hours, minutes, seconds = match.groups()
    total = 0.0
    if hours:
        total += int(hours) * 3600
    if minutes:
        total += int(minutes) * 60
    if seconds:
        total += float(seconds)
    return total


def parse_size(value):
    """Convert Nextflow byte-size strings (e.g. '5.4 GB') to gigabytes."""
    if pd.isna(value) or value in ("-", ""):
        return np.nan
    num, unit = value.strip().split()
    return float(num) * _SIZE_UNITS[unit] / 1e9


def parse_pct(value):
    """Convert a '%cpu' string (e.g. '198.0%') to a float."""
    if pd.isna(value) or value in ("-", ""):
        return np.nan
    return float(value.strip().rstrip("%"))


# -- load + tidy ------------------------------------------------------------

def load_trace(path):
    df = pd.read_csv(path, sep="\t")
    df["stage"] = df["name"].map(broad_stage)
    df = df.dropna(subset=["stage"]).copy()
    df["process"] = df["name"].map(extract_process)

    df["realtime_min"] = df["realtime"].map(parse_duration) / 60
    df["cpu_pct"] = df["%cpu"].map(parse_pct)
    df["peak_rss_gb"] = df["peak_rss"].map(parse_size)
    df["io_gb"] = df["rchar"].map(parse_size) + df["wchar"].map(parse_size)
    return df


def process_sequence():
    """Full ordered (stage, process) list, pipeline order, for x-axis/table ordering."""
    return [(stage, proc) for stage in STAGE_ORDER for proc in PROCESS_ORDER.get(stage, [])]


# -- plotting ---------------------------------------------------------------

METRICS = [
    ("realtime_min", "Wall-clock time (min)"),
    ("cpu_pct", "CPU utilisation (%)"),
    ("peak_rss_gb", "Peak memory, RSS (GB)"),
    ("io_gb", "Total I/O, read + write (GB)"),
]

STAGE_BASE_COLOR = {
    "QUALITY_CONTROL": "#4C72B0",
    "REFERENCE_ALIGNMENT": "#DD8452",
    "VARIANT_CALLING": "#55A868",
    "JOINT_GENOTYPING": "#C44E52",
    "VARIANT_FILTERING": "#8172B2",
    "VARIANT_EFFECT_PREDICTION": "#64B5CD",
    "MULTIQC": "#DA8BC3",
}

MIN_TASKS_FOR_BOX = 2  # process groups below this (e.g. single-task MultiQC) are dropped

PANEL_LETTERS = "ABCD"


def _shades(base_hex, n):
    """n distinct lightness shades of base_hex (same hue/saturation), for processes within a stage."""
    r, g, b = mcolors.to_rgb(base_hex)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if n == 1:
        lightness_values = [l]
    else:
        lightness_values = np.linspace(max(l - 0.28, 0.20), min(l + 0.28, 0.82), n)
    return [colorsys.hls_to_rgb(h, lig, s) for lig in lightness_values]


def _fmt_value(y):
    """Plain-decimal formatting (no scientific notation) sized to the value's magnitude."""
    if y >= 100:
        return f"{y:.0f}"
    if y >= 10:
        return f"{y:.1f}"
    if y >= 1:
        return f"{y:.2f}"
    return f"{y:.3f}"


def _annotate_medians(ax, bp):
    """Label each box with its median value, anchored a fixed screen-distance above its
    own upper whisker (constant point-offset, not scale-dependent), so labels stay
    attached to their box regardless of how small or large that box is relative to others."""
    ylim = ax.get_ylim()
    span = ylim[1] - ylim[0]
    ax.set_ylim(ylim[0], ylim[1] + 0.12 * span)

    for i, med_line in enumerate(bp["medians"]):
        median = med_line.get_ydata()[0]
        cap_top = bp["caps"][2 * i + 1].get_ydata()[0]  # upper whisker cap
        ax.annotate(
            _fmt_value(median),
            xy=(i + 1, cap_top), xycoords="data",
            xytext=(0, 4), textcoords="offset points",
            ha="center", va="bottom", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8),
        )


def plot_metrics(df, out_path):
    counts = df.groupby(["stage", "process"]).size()
    groups = [(s, p) for s, p in process_sequence() if counts.get((s, p), 0) >= MIN_TASKS_FOR_BOX]
    labels = [p.replace("__", "\n") for _, p in groups]

    colors = []
    for stage in STAGE_ORDER:
        procs_in_stage = [p for s, p in groups if s == stage]
        colors.extend(_shades(STAGE_BASE_COLOR[stage], len(procs_in_stage)))

    fig, axes = plt.subplots(4, 1, figsize=(16, 18))

    for ax, letter, (col, ylabel) in zip(axes, PANEL_LETTERS, METRICS):
        data = [df.loc[(df["stage"] == s) & (df["process"] == p), col].dropna().values for s, p in groups]
        bp = ax.boxplot(
            data,
            tick_labels=labels,
            showfliers=True,
            flierprops=dict(marker="o", markersize=3, alpha=0.4, markeredgecolor="none"),
            medianprops=dict(color="black"),
            patch_artist=True,
            widths=0.6,
        )
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.85)
        for flier, color in zip(bp["fliers"], colors):
            flier.set_markerfacecolor(color)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.tick_params(axis="y", labelsize=10.5)
        ax.grid(axis="y", linestyle="--", alpha=0.3)
        _annotate_medians(ax, bp)
        ax.text(-0.045, 1.03, letter, transform=ax.transAxes, fontsize=17, fontweight="bold")

        # Shade stage groupings consistently across all panels.
        start = 0
        for stage in STAGE_ORDER:
            n = sum(1 for s, _ in groups if s == stage)
            if n:
                ax.axvspan(start + 0.5, start + n + 0.5, color=STAGE_BASE_COLOR[stage], alpha=0.06, zorder=0)
                start += n

    # Stage name headers above the top panel.
    start = 0
    for stage in STAGE_ORDER:
        n = sum(1 for s, _ in groups if s == stage)
        if n:
            axes[0].text(
                start + (n + 1) / 2, 1.14, STAGE_LABELS[stage],
                transform=axes[0].get_xaxis_transform(),
                ha="center", va="bottom", fontsize=12.5, fontweight="bold",
                color=STAGE_BASE_COLOR[stage],
            )
            start += n

    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right", fontsize=12)

    fig.suptitle("Nextflow task resource usage by process, grouped by pipeline stage", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


def _summarize(series):
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    return pd.Series({
        "count": series.count(),
        "median": series.median(),
        "mean": series.mean(),
        "std": series.std(),
        "min": series.min(),
        "max": series.max(),
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
    })


def write_stats(df, out_path):
    metrics = [col for col, _ in METRICS]
    frames = {
        col: df.groupby(["stage", "process"])[col].apply(_summarize).unstack()
        for col in metrics
    }
    summary = pd.concat(frames, axis=1)

    order = process_sequence()
    order = [t for t in order if t in summary.index]
    remaining = [t for t in summary.index if t not in order]
    summary = summary.reindex(order + sorted(remaining))
    summary.to_csv(out_path, sep="\t")
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="Nextflow execution_trace TSV file")
    parser.add_argument("-o", "--output", default="execution_trace_metrics.png", help="Output figure path")
    parser.add_argument("--stats", default=None, help="Optional output TSV of per-stage summary statistics")
    args = parser.parse_args()

    df = load_trace(args.input)
    plot_metrics(df, args.output)
    if args.stats:
        write_stats(df, args.stats)


if __name__ == "__main__":
    main()
