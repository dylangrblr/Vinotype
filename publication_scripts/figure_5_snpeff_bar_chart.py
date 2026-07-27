#!/usr/bin/env python3
"""
Plot SnpEff annotation distribution as a grouped bar chart comparing
Benchmark_shared and Vinotype_private datasets.

Raw SnpEff region types are collapsed into six biologically meaningful
groups ordered from most to least genic. Percentages are recalculated
relative to total annotation events within each dataset (not unique
variants, since one variant can carry multiple annotations). Counts 
are derived from SnpEff html output

Usage:
    figure_5_snpeff_bar_chart.py -o test
"""

# ==================== IMPORTS ================================================#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import argparse

# ==================== RAW DATA ===============================================#

# Raw counts as reported by SnpEff (annotation events, not unique variants).
# One variant can receive multiple annotations.
RAW_COUNTS = {
    "Benchmark_shared": {
        "DOWNSTREAM":               6_313,
        "EXON":                     3_157,
        "INTERGENIC":               1_827,
        "INTRON":                   4_613,
        "SPLICE_SITE_ACCEPTOR":         1,
        "SPLICE_SITE_DONOR":            2,
        "SPLICE_SITE_REGION":         202,
        "TRANSCRIPT":                  76,
        "UPSTREAM":                 5_436,
        "UTR_3_PRIME":                742,
        "UTR_5_PRIME":                616,
    },
    "Vinotype_private": {
        "DOWNSTREAM":           9_119_986,
        "EXON":                 1_509_156,
        "INTERGENIC":           9_663_542,
        "INTRON":               4_843_677,
        "SPLICE_SITE_ACCEPTOR":     3_678,
        "SPLICE_SITE_DONOR":        2_865,
        "SPLICE_SITE_REGION":      63_115,
        "TRANSCRIPT":             230_698,
        "UPSTREAM":             9_323_161,
        "UTR_3_PRIME":            743_965,
        "UTR_5_PRIME":            511_157,
    },
}

DATASET_ORDER = ["Benchmark_shared", "Vinotype_private"]

# Ts/Tv ratios reported by SnpEff for each dataset
TSTV_RATIOS = {
    "Benchmark_shared": 3.569,
    "Vinotype_private":  2.320,
}

# Visual style per dataset: solid fill for Benchmark, hatched for Vinotype
DATASET_HATCHES = {
    "Benchmark_shared": "",
    "Vinotype_private":  "////",
}

# Group definitions:
#   Exonic regions → all transcribed sequence (exons, UTRs, intragenic)
#   Intronic       → intronic variants
#   Splice sites   → all three splice site categories
#   Upstream       → within 5 kb upstream of TSS
#   Downstream     → within 5 kb downstream of TES
#   Intergenic     → distal from any annotated gene
GROUP_DEFINITIONS = {
    "Exonic\nregions": [
        "EXON",
        "TRANSCRIPT",
        "UTR_3_PRIME",
        "UTR_5_PRIME",
    ],
    "Intronic": [
        "INTRON",
    ],
    "Splice\nsites": [
        "SPLICE_SITE_ACCEPTOR",
        "SPLICE_SITE_DONOR",
        "SPLICE_SITE_REGION",
    ],
    "Upstream": [
        "UPSTREAM",
    ],
    "Downstream": [
        "DOWNSTREAM",
    ],
    "Intergenic": [
        "INTERGENIC",
    ],
}

# Display order: most genic (left) --> most distal (right)
GROUP_ORDER = [
    "Exonic\nregions",
    "Splice\nsites",
    "Intronic",
    "Upstream",
    "Downstream",
    "Intergenic",
]

# Gradient from warm/saturated (genic) to cool/grey (distal)
#GROUP_COLORS = {
#    "Exonic\nregions": "#2e7d32",  # dark green  – most functionally relevant
#    "Splice\nsites":   "#f57f17",  # amber       – high impact, small count
#    "Intronic":        "#81c784",  # light green – genic but non-coding
#    "Upstream":        "#90caf9",  # light blue  – potential regulatory role
#    "Downstream":      "#64b5f6",  # mid blue    – potential regulatory role
#    "Intergenic":      "#b0bec5",  # grey        – most distal
#}

GROUP_COLORS = {
    "Exonic\nregions": "#8d2f2f",  # dark red  – most functionally relevant
    "Splice\nsites":   "#ec5552",  # medium red       – high impact, small count
    "Intronic":        "#ee8590",  # light red – genic but non-coding
    "Upstream":        "#2e7d32",  # dark green  – potential regulatory role
    "Downstream":      "#a5d6a7",  # light green – potential regulatory role
    "Intergenic":      "#b0bec5",  # grey        – most distal
}


# ==================== HELPER FUNCTIONS =======================================#


def build_groups(raw_counts, group_definitions):
    """
    Aggregate raw SnpEff counts into simplified biological groups.

    Parameters
    ----------
    raw_counts        : dict  {category_name: int}
    group_definitions : dict  {group_label: [list of raw category names]}

    Returns
    -------
    dict  {group_label: aggregated_count}
    """
    grouped = {}
    for group_label, categories in group_definitions.items():
        # Sum counts for all raw categories belonging to this group.
        # .get() returns 0 for any category not found in raw_counts
        grouped[group_label] = sum(
            raw_counts.get(cat, 0) for cat in categories
        )
    return grouped


def compute_pct(grouped_counts):
    """
    Recalculate percentages from aggregated group counts.

    Percentages are relative to the total annotation events across all
    groups, NOT the number of unique variants.

    Parameters
    ----------
    grouped_counts : dict  {group_label: int}

    Returns
    -------
    dict  {group_label: float (percentage)}
    """
    total = sum(grouped_counts.values())
    return {
        label: (count / total * 100)
        for label, count in grouped_counts.items()
    }


# ==================== MAIN ===================================================#


def main(output_prefix):

    # ------------------------------------------------------------------
    # 1. Aggregate raw counts into groups and recalculate percentages
    #    for each dataset
    # ------------------------------------------------------------------
    all_grouped_counts = {}
    all_grouped_pct = {}
    for dataset in DATASET_ORDER:
        all_grouped_counts[dataset] = build_groups(RAW_COUNTS[dataset], GROUP_DEFINITIONS)
        all_grouped_pct[dataset] = compute_pct(all_grouped_counts[dataset])

    # ------------------------------------------------------------------
    # 2. Export grouped summary to TSV
    # ------------------------------------------------------------------
    rows = []
    for dataset in DATASET_ORDER:
        for group in GROUP_ORDER:
            rows.append({
                "dataset": dataset,
                "group":   group.replace("\n", " "),
                "count":   all_grouped_counts[dataset][group],
                "percent": all_grouped_pct[dataset][group],
            })
    df = pd.DataFrame(rows)
    tsv_out = f"{output_prefix}_snpeff_grouped_summary.tsv"
    df.to_csv(tsv_out, sep="\t", index=False)
    print(f"📄 Summary table saved: {tsv_out}")

    # ------------------------------------------------------------------
    # 3. Plot grouped bars
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 6))

    n_groups = len(GROUP_ORDER)
    bar_width = 0.35
    x = np.arange(n_groups)

    # Pre-compute  axis ceiling and label gap before plot loop
    max_pct      = max(
        pct
        for ds in DATASET_ORDER
        for pct in all_grouped_pct[ds].values()
    )
    label_offset = max_pct * 0.01   # 3 % of the tallest bar

    for i, dataset in enumerate(DATASET_ORDER):
        offset = (i - 0.4) * bar_width          # centre around each tick
        pcts   = [all_grouped_pct[dataset][g] for g in GROUP_ORDER]
        colors = [GROUP_COLORS[g] for g in GROUP_ORDER]
        hatch  = DATASET_HATCHES[dataset]

        bars = ax.bar(
            x + offset,
            pcts,
            bar_width,
            color=colors,
            hatch=hatch,
            edgecolor="black",
            linewidth=0.8,
        )

        # Annotate each bar with its rounded percentage
        for bar, pct in zip(bars, pcts):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + label_offset,
                f"{pct:.2f}%",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )

    # ------------------------------------------------------------------
    # 4. Legend: dataset patterns first, then group colours
    # ------------------------------------------------------------------
    dataset_handles = [
        Patch(
            facecolor="white",
            edgecolor="black",
            hatch=DATASET_HATCHES[ds],
            label=ds,
        )
        for ds in DATASET_ORDER
    ]
    group_handles = [
        Patch(
            facecolor=GROUP_COLORS[g],
            edgecolor="black",
            label=g.replace("\n", " "),
        )
        for g in GROUP_ORDER
    ]
    # Blank separator between the two legend sections
    separator = Patch(facecolor="none", edgecolor="none", label="")

    ax.legend(
        handles=dataset_handles + [separator] + group_handles,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        borderaxespad=0,
        fontsize=11,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(GROUP_ORDER, fontsize=13)
    ax.set_ylabel("Percentage of annotation events (%)", fontsize=13)
    ax.set_xlabel("Genomic region", fontsize=13)
    ax.set_title(
        "Distribution of SnpEff Annotations by Genomic Region",
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylim(0, max_pct * 1.18)
    ax.margins(x=0.05)

    # Footer note
    fig.text(
        0.5,
        -0.02,
        "Note: bars represent % of total annotation events per dataset; "
        "one variant may carry multiple annotations.",
        ha="center",
        fontsize=11,
        color="grey",
        style="italic",
    )

    plt.tight_layout()

    # ------------------------------------------------------------------
    # 5. Ts/Tv ratio
    # ------------------------------------------------------------------
    fig.canvas.draw()
    renderer        = fig.canvas.get_renderer()
    legend          = ax.get_legend()
    legend_bbox     = legend.get_window_extent(renderer=renderer)
    ax_bbox         = ax.get_window_extent(renderer=renderer)
    legend_bottom_y = (legend_bbox.y0 - ax_bbox.y0) / ax_bbox.height

    legend_frame = legend.get_frame()
    frame_fc = legend_frame.get_facecolor()
    frame_ec = legend_frame.get_edgecolor()
    frame_lw = legend_frame.get_linewidth()

    tstv_lines = ["Ts/Tv ratio"] + [
        f"  {ds}: {TSTV_RATIOS[ds]:.2f}" for ds in DATASET_ORDER
    ]
    ax.text(
        1.01,
        legend_bottom_y - 0.03,
        "\n".join(tstv_lines),
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=11,
        fontfamily="monospace",
        bbox=dict(boxstyle="square,pad=0.4", facecolor=frame_fc, edgecolor=frame_ec, linewidth=frame_lw),
    )

    png_out = f"{output_prefix}_snpeff_distribution.png"
    plt.savefig(png_out, dpi=800, bbox_inches="tight")
    print(f"✅ Plot saved: {png_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Plot SnpEff annotation distribution collapsed into "
            "six genomic region groups."
        )
    )
    parser.add_argument(
        "-o",
        "--out",
        default="output",
        help="Output file prefix [default: output]",
    )
    args = parser.parse_args()
    main(args.out)
