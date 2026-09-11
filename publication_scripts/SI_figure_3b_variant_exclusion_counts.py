#!/usr/bin/env python3
"""Plot variant exclusion counts by filter or annotation."""

from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FILTERS = [
    ("Strand odds ratio (SOR)", "≥ 3.0", 8),
    ("Fisher strand bias (FS)", "≥ 60.0", 5),
    ("Read position rank sum (ReadPosRankSum)", "≤ −8.0", 0),
    ("Mapping quality rank sum (MQRankSum)", "≤ −12.5", 0),
    ("Quality by depth (QD)", "≤ 2.0", 2),
    ("Allele count – singletons (AC)", "< 2", 6),
    ("Read depth (DP)", "> 10% of samples with DP < 5", 493),
    ("Genotype quality (GQ)", "> 25% of samples with GQ ≤ 10", 99),
    ("Allele depth (AD) – AR", "> 10% of AR have AD < 2 reads", 34),
    ("Allele depth (AD) – AA", "> 10% of AA have AD < 4 reads", 546),
    ("Variant allele frequency (VAF) – AR", "> 80% of AR have a VAF of < 0.5", 33),
    ("Genotype call rate", "< 80% completeness", 91),
]

COLOR = "#4C72B0"


def plot_exclusion_counts(output_path: str) -> None:
    """Create and save the horizontal exclusion-count barplot."""
    filters = [f"{item[0]}\n[{item[1]}]" for item in FILTERS]
    counts = np.array([item[2] for item in FILTERS])
    positions = np.arange(len(filters))

    fig, ax = plt.subplots(figsize=(14, 8.5))
    bars = ax.barh(
        positions,
        counts,
        color=COLOR,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.85,
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(filters, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Number of excluded variants", fontsize=11, fontweight="bold")
    ax.set_title(
        "Variant exclusion counts by filter or annotation",
        fontsize=13,
        fontweight="bold",
        pad=14,
    )
    ax.text(
        0.0,
        1.01,
        "SNP variant class",
        transform=ax.transAxes,
        fontsize=10,
        color="dimgray",
        style="italic",
    )

    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(True, alpha=0.3, axis="x")
    ax.set_axisbelow(True)
    ax.set_xlim(0, counts.max() * 1.18)

    for bar, count in zip(bars, counts):
        y_center = bar.get_y() + bar.get_height() / 2
        ax.text(
            count + counts.max() * 0.012,
            y_center,
            f"{count:,}",
            ha="left",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=COLOR,
        )

    ax.text(
        0.0,
        -0.08,
        "Thresholds are shown in brackets; values are the supplied exclusion counts.",
        transform=ax.transAxes,
        fontsize=9,
        color="dimgray",
        style="italic",
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Plot saved → {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot SNP variant exclusion counts by filter or annotation."
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variant_exclusion_counts.png",
        help="Output PNG path (default: variant_exclusion_counts.png)",
    )
    args = parser.parse_args()
    plot_exclusion_counts(args.output)


if __name__ == "__main__":
    main()