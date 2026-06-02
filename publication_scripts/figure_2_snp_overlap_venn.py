#!/usr/bin/env python3
"""
Venn diagram showing SNP overlap between three datasets:
  - SNP-array benchmark  (10,228 SNPs)
  - Vinotype raw          (41,943,184 SNPs)
  - Vinotype filtered     (15,846,099 SNPs)

All datasets cover the same 220 samples.

Assumptions:
  - Filtered is a strict subset of Raw (filtered SNPs ⊆ raw SNPs).
  - Array ∩ Raw  = 9,507   (array SNPs recovered by raw calls)
  - Array ∩ Filtered = 9,434   (array SNPs recovered by filtered calls)
  - Therefore Array ∩ Raw ∩ Filtered = 9,434
  - Array ∩ Raw only (not Filtered) = 9,507 − 9,434 = 73
  - Array only = 10,228 − 9,507 = 721
  - Filtered only = 0  (Filtered ⊆ Raw, so no Filtered-exclusive SNPs outside Raw)
  - Raw ∩ Filtered (not Array) = 15,846,099 − 9,434 = 15,836,665
  - Raw only (not Array, not Filtered) = 41,943,184 − 9,507 − 15,836,665 = 26,097,012
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Patch
from matplotlib.colors import to_rgba
import numpy as np

# ── helpers ──────────────────────────────────────────────────────────────────

def fmt(n):
    """Format large numbers with comma separators."""
    return f"{n:,}"

# ── counts ───────────────────────────────────────────────────────────────────

array_total     = 10_228
raw_total       = 41_943_184
filtered_total  = 15_846_099

array_raw_filt  = 9_434          # all three
array_raw_only  = 9_507 - 9_434  # array ∩ raw, NOT filtered  = 73
array_only      = array_total - 9_507  # 721

raw_filt_only   = filtered_total - array_raw_filt   # raw ∩ filtered, NOT array = 15,836,665
raw_only        = raw_total - 9_507 - raw_filt_only  # raw only = 26,097,012
# filtered_only = 0  (Filtered ⊆ Raw)
# array_filt_only = 0  (Filtered ⊆ Raw)

# ── figure layout ────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(11, 8))
ax.set_xlim(0, 11)
ax.set_ylim(0, 8)
ax.set_aspect("equal")
ax.axis("off")

# Ellipse parameters  (cx, cy, width, height, angle)
# Array   – top centre, smaller circle
# Raw     – bottom left, large
# Filtered– bottom right, medium

alpha = 0.28

circles = {
    "Array":    dict(xy=(5.5, 5.2), width=3.2, height=2.8, angle=0,
                     color="#F4C97A", label="SNP-Array\nBenchmark"),
    "Raw":      dict(xy=(3.8, 3.4), width=5.6, height=4.0, angle=-10,
                     color="#89B4D9", label="Vinotype\nRaw"),
    "Filtered": dict(xy=(7.0, 3.4), width=4.2, height=3.4, angle=10,
                     color="#85C9A0", label="Vinotype\nFiltered"),
}

for key, p in circles.items():
    e = Ellipse(xy=p["xy"], width=p["width"], height=p["height"],
                angle=p["angle"],
                facecolor=to_rgba(p["color"], alpha),
                edgecolor=p["color"], linewidth=1.2, zorder=2)
    ax.add_patch(e)

# ── annotation positions & values ────────────────────────────────────────────

annotations = [
    # (x,    y,    text,                                          fontsize, ha)
    (5.5,  6.55,  f"Array only\n{fmt(array_only)}",               9,  "center"),
    (2.2,  2.1,   f"Raw only\n{fmt(raw_only)}",                   9,  "center"),
    (9.1,  2.1,   f"Filtered only\n0",                            9,  "center"),
    (4.0,  5.05,  f"Array ∩ Raw\n(not Filtered)\n{fmt(array_raw_only)}", 8.5, "center"),
    (5.5,  4.35,  f"Array ∩ Raw\n∩ Filtered\n{fmt(array_raw_filt)}", 8.5, "center"),
    (6.5,  2.25,  f"Raw ∩ Filtered\n(not Array)\n{fmt(raw_filt_only)}", 8.5, "center"),
    # Array ∩ Filtered only = 0, omit for clarity
]

for x, y, text, fs, ha in annotations:
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs,
            fontweight="bold", color="#1a1a1a", zorder=5,
            linespacing=1.4)

# ── legend labels ─────────────────────────────────────────────────────────────

legend_entries = [
    Patch(facecolor=to_rgba("#F4C97A", 0.55), edgecolor="#F4C97A",
          linewidth=1.2, label=f"SNP-Array Benchmark\n({fmt(array_total)} SNPs)"),
    Patch(facecolor=to_rgba("#89B4D9", 0.55), edgecolor="#89B4D9",
          linewidth=1.2, label=f"Vinotype Raw\n({fmt(raw_total)} SNPs)"),
    Patch(facecolor=to_rgba("#85C9A0", 0.55), edgecolor="#85C9A0",
          linewidth=1.2, label=f"Vinotype Filtered\n({fmt(filtered_total)} SNPs)"),
]

ax.legend(handles=legend_entries, loc="lower center",
          bbox_to_anchor=(0.5, -0.02), ncol=3,
          fontsize=9, framealpha=0.85, edgecolor="#cccccc",
          handlelength=1.8, handleheight=1.4)

# ── title ─────────────────────────────────────────────────────────────────────

ax.set_title("Overlap Among Vinotype-Predicted and Benchmark SNPs",
             fontsize=14, fontweight="bold", pad=12, color="#1a1a1a")

plt.tight_layout()
out_path = "docs/snp_overlap_venn.png"
fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
print(f"Saved → {out_path}")
plt.close(fig)
