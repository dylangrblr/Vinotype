#!/usr/bin/env python3
"""
Remove redundant GO terms from gProfiler results using goatools.

Requires gprofiler output tsv (final_go_gprofiler.tsv) and the GO DAG (go-basic.obo).
-------
"""

import subprocess

import pandas as pd
from goatools.obo_parser import GODag


# ── Parsing ──────────────────────────────────────────────────────────────────

def parse_gprofiler_tsv(filepath: str) -> pd.DataFrame:
    """Parse a gProfiler result TSV into a DataFrame, retaining only GO terms."""
    df = pd.read_csv(filepath, sep="\t")
    df.columns = df.columns.str.strip()
    df = df[df["source"].str.startswith("GO:")].copy()
    return df


# ── Redundancy detection ──────────────────────────────────────────────────────

def get_all_ancestors(go_id: str, godag: GODag) -> set[str]:
    """Collect all ancestor GO IDs by traversing up the DAG."""
    ancestors = set()
    to_visit = list(godag[go_id].parents) if go_id in godag else []
    while to_visit:
        parent = to_visit.pop()
        if parent.item_id not in ancestors:
            ancestors.add(parent.item_id)
            to_visit.extend(parent.parents)
    return ancestors


def find_redundant_ancestors(
    go_terms: list[str], godag: GODag
) -> tuple[set[str], set[str]]:
    """
    Identify terms that are ancestors of other terms in the same list.
    Returns (terms_to_remove, terms_to_keep).
    """
    go_set = set(go_terms)
    valid_terms = {term for term in go_set if term in godag}
    missing_terms = go_set - valid_terms
    terms_to_remove = set()

    for term in sorted(missing_terms):
        print(f"  Warning: {term} not found in OBO, skipping.")

    for term in valid_terms:
        ancestors = get_all_ancestors(term, godag)
        terms_to_remove.update(ancestors & valid_terms)

    terms_to_keep = valid_terms - terms_to_remove
    print(
        f"  Redundant ancestors: {len(terms_to_remove)} | "
        f"Retained: {len(terms_to_keep)} | "
        f"Missing from OBO: {len(missing_terms)}"
    )
    return terms_to_remove, terms_to_keep


def get_redundant_pairs(
    terms_to_remove: set[str],
    terms_to_keep: set[str],
    godag: GODag,
) -> set[str]:
    """
    Return only the subset of terms directly involved in a
    parent-child redundancy within the user's list.
    Useful for a focused, readable plot.
    """
    pairs: set[str] = set()
    for term in terms_to_keep:
        if term not in godag:
            continue
        ancestors = get_all_ancestors(term, godag)
        shared = ancestors & terms_to_remove
        if shared:
            pairs.add(term)
            pairs.update(shared)
    return pairs


# ── Plotting ──────────────────────────────────────────────────────────────────

# Lightbright palette colors (http://www.color-hex.com/color-palette/51895)
COLOR_RETAINED  = "#e6fad2"  # light green  — most specific term, retained
COLOR_REDUNDANT = "#fad2d2"  # light pink   — ancestor of another term, removed


def write_go_color_file(
    go_terms: list[str],
    terms_to_remove: set[str],
    terms_to_keep: set[str],
    godag: GODag,
    filepath: str,
) -> None:
    """Write a goatools go_plot color file distinguishing retained vs redundant terms."""
    lines = ["# goatools go_plot color file — retained=green, redundant=pink"]
    for term in sorted(go_terms):
        if term not in godag:
            continue
        name = godag[term].name
        if term in terms_to_remove:
            lines.append(f"{COLOR_REDUNDANT} {term}  # {name} (redundant)")
        elif term in terms_to_keep:
            lines.append(f"{COLOR_RETAINED} {term}  # {name} (retained)")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Color file written: {filepath}")


def plot_go_hierarchies(
    color_file: str,
    output_file: str = "go_hierarchy.png",
    obo_file: str = "go-basic.obo",
) -> None:
    """Call goatools go_plot CLI with the generated color file."""
    cmd = [
        "goatools", "go_plot",
        f"--go_file={color_file}",
        f"--obo={obo_file}",
        "-o", output_file,
    ]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  go_plot stderr:\n{result.stderr}")
    else:
        print(f"  Saved: {output_file}")


# ── Main ──────────────────────────────────────────────────────────────────────

# 1. Load GO DAG
print("Loading GO DAG...")
godag = GODag("go-basic.obo")

# 2. Parse gProfiler output
df = parse_gprofiler_tsv("final_go_gprofiler.tsv")
print(f"Significant GO terms loaded from gProfiler: {len(df)}")
print(df.groupby("source").size().to_string())

go_terms = df["native"].tolist()

# 3. Optionally pre-filter by source before redundancy removal, e.g.:
# df = df[df["source"] == "GO:BP"]
# go_terms = df["native"].tolist()
# print(f"After source filter: {len(go_terms)} terms")

# 4. Find redundant ancestors
print("\nIdentifying redundant ancestors...")
terms_to_remove, terms_to_keep = find_redundant_ancestors(go_terms, godag)

# 5. Write color file and plot
print("\nGenerating hierarchy plot...")
write_go_color_file(
    go_terms, terms_to_remove, terms_to_keep, godag,
    filepath="gprofiler_go_colors.txt",
)
plot_go_hierarchies(
    color_file="gprofiler_go_colors.txt",
    output_file="gprofiler_go_before_filter.png",
)

# 6. Filter DataFrame and save
df_final = df[df["native"].isin(terms_to_keep)].copy()
print(f"\nFinal dataset: {len(df_final)} terms retained")
print(df_final[["source", "native", "name", "p_value"]].to_string(index=False))

df_final.to_csv("gprofiler_redundancy_filtered.tsv", sep="\t", index=False)
print("\nSaved: gprofiler_redundancy_filtered.tsv")