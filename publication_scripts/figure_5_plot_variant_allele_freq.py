#!/usr/bin/env python3
"""
Plot per-variant allele frequencies for a given gene from a SnpEff annotated VCF,
and annotate samples carrying a selected amino acid substitution.

Usage:
   python figure_5_plot_variant_allele_freq.py -v file.vcf.gz -g AGL11 -m Arg197Leu
"""

# -- imports --

from cyvcf2 import VCF
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import argparse
import re
import os

# -- helpers --

def parse_ann_field(ann_field):
    """Extract amino acid change and predicted impact from ANN string."""

    #If annotation field is empty, return None, None
    if ann_field is None:
        return None, None

    annotations = ann_field.split(",")
    for ann in annotations:
        fields = ann.split("|")
        aa_field = next((f for f in fields if f.startswith("p.")), None)
        impact_field = next(
            (f for f in fields if f in ["HIGH", "MODERATE", "LOW", "MODIFIER"]),
            None,
        )
        if aa_field:
            aa_change = re.sub(r"p\\.?", "", aa_field)
            return aa_change, impact_field
    return None, impact_field

def compute_allele_frequency(genotypes):
    """Compute simple allele frequency (het = 0.5, hom-alt = 1)."""
    het = sum(gt == 1 for gt in genotypes)
    hom_alt = sum(gt == 3 for gt in genotypes)
    total = len(genotypes)
    return (het * 0.5 + hom_alt) / total if total > 0 else None


def normalize_mutation_name(mutation_name):
    """Strip an optional 'p.' prefix so user input matches ANN amino-acid calls."""
    if mutation_name is None:
        return None
    return re.sub(r"^p\.?", "", mutation_name)

# -- main --

def main(vcf_file, gene_name, mutation_name, output_prefix):

    vcf = VCF(vcf_file)
    samples = vcf.samples
    records = []
    carriers = []
    mutation_query = mutation_name

    for variant in vcf:
        ann = variant.INFO.get("ANN")
        if ann and gene_name in ann:
            aa_change, impact = parse_ann_field(ann)
            if not aa_change:
                continue
            af = compute_allele_frequency(variant.gt_types)
            if af is None:
                continue
            records.append(
                {
                    "pos": variant.POS,
                    "aa_change": aa_change,
                    "Impact": impact,
                    "allele_freq": af,
                }
            )
            if mutation_query and mutation_query in aa_change:
                for sample, gt in zip(samples, variant.gt_types):
                    if gt in (1, 3):  # het or hom-alt
                        carriers.append(sample)

    if not records:
        print(f"No variants found for {gene_name}")
        return

    df = pd.DataFrame(records).sort_values("allele_freq", ascending=False).reset_index(drop=True)

# -- plot --

    plt.figure(figsize=(14, 6))
    ax = plt.gca()

    impact_palette = {
        "HIGH": "#e64a19",
        "MODERATE": "#234e24",
        "LOW": "#81c784",
        "MODIFIER": "#b0bec5",
    }

    colors = df["Impact"].map(impact_palette).fillna("#999999")
    x_positions = list(range(len(df)))
    bars = ax.bar(
        x_positions,
        df["allele_freq"],
        color=colors,
        width=0.8,
        edgecolor="black",
        linewidth=0.8,
    )

    legend_handles = [
        Patch(facecolor=impact_palette[label], edgecolor="black", label=label)
        for label in ["MODERATE", "LOW"]
    ]
    ax.legend(handles=legend_handles, title="Impact")

    # Highlight the selected mutation bar with a distinct, colorblind-friendly color.
    n_carriers = len(sorted(set(carriers)))
    target_rows = df.index[df["aa_change"] == mutation_query].tolist()
    if target_rows:
        target_idx = target_rows[0]
        if target_idx < len(bars):
            target_patch = bars[target_idx]
            target_patch.set_facecolor("#d7301f")
            target_patch.set_edgecolor("black")
            target_patch.set_linewidth(1.2)

            bar_center_x = target_patch.get_x() + (target_patch.get_width() / 2)
            bar_top_y = target_patch.get_height()
            ax.annotate(
                f"{mutation_query}\nn = {n_carriers} carriers",
                xy=(bar_center_x, bar_top_y),
                xytext=(0, 18),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                color="#d7301f",
                arrowprops={"arrowstyle": "->", "color": "#d7301f", "lw": 1.3},
                bbox={"boxstyle": "round,pad=0.25", "fc": "#fff5f0", "ec": "#d7301f", "alpha": 0.95},
                annotation_clip=False,
            )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(df["aa_change"], rotation=45, ha="right", fontsize=8)
    ax.tick_params(axis="x", pad=8)
    ax.set_xlim(-0.55, len(df) - 0.45)
    plt.xlabel("Amino Acid Substitution")
    plt.ylabel("Allele Frequency")
    plt.title(f"Variant Spectrum in {gene_name}")
    ax.margins(y=0.12)

    # Export carrier names to a supplementary table; keep the figure uncluttered.
    if carriers:
        unique_carriers = sorted(set(carriers))

        carriers_df = pd.DataFrame({"sample": unique_carriers})
        carriers_out = f"{output_prefix}_{gene_name}_{mutation_query}_carriers.tsv"
        carriers_df.to_csv(carriers_out, sep="\t", index=False)

        print(f"🧬 Carriers of {mutation_query}: {', '.join(unique_carriers)}")
        print(f"🗂️ Carrier table saved: {carriers_out}")
    else:
        print(f"No carriers of {mutation_query} found in {gene_name} within {vcf_file}")

    plt.tight_layout()
    out_png = f"{output_prefix}_{gene_name}_variants_{mutation_query}.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"✅ Plot saved: {out_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot variant frequencies for a gene and label carriers of a specific amino acid substitution."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input annotated VCF file")
    parser.add_argument("-g", "--gene", required=True, help="Gene name (e.g., AGL11)")
    parser.add_argument(
        "-m", "--mutation", required=True, help="Amino acid change, e.g. Arg197Leu"
    )
    parser.add_argument(
        "-o",
        "--out",
        default="output",
        help="Output file prefix [default: output]",
    )
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"VCF not found: {args.vcf}")

    main(args.vcf, args.gene, args.mutation, args.out)