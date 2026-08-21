#!/usr/bin/env python3
"""
ESM2 Variant Effect Predictor — HPC/container edition.

Usage:
    python esm_vep.py \
        --fasta_dir /data/fasta \
        --var_dir   /data/variants \
        --output    /data/results/scores.csv
"""
import os
import re
import sys
import argparse
import warnings
import logging

import torch
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm
from transformers import AutoTokenizer, EsmForMaskedLM
from transformers import logging as hf_logging

warnings.filterwarnings("ignore")
hf_logging.set_verbosity_error()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def find_matching_files(fasta_dir: str, var_dir: str) -> list[tuple]:
    matches = []
    for fname in os.listdir(fasta_dir):
        if not fname.endswith(".fasta"):
            continue
        base = fname.replace(".fasta", "")
        var_path = os.path.join(var_dir, f"{base}.var.txt")
        if os.path.exists(var_path):
            matches.append((fname, var_path, base))
        else:
            logger.warning("No variant file for %s — skipping", base)
    return matches


def load_protein_sequence(path: str) -> str | None:
    for rec in SeqIO.parse(path, "fasta"):
        return str(rec.seq)
    return None


def load_variants(path: str) -> list[tuple]:
    variants = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"^([A-Z])(\d+)([A-Z])$", line)
            if m:
                variants.append((m.group(1), int(m.group(2)), m.group(3)))
    return variants


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def load_model(model_name: str, device: torch.device):
    logger.info("Loading model: %s", model_name)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name, ignore_mismatched_sizes=True)
    return model.to(device).eval(), tokenizer


def compute_variant_score(
    sequence: str,
    variant: tuple,
    model: EsmForMaskedLM,
    tokenizer: AutoTokenizer,
    device: torch.device,
) -> float:
    wt, pos, mt = variant
    idx0 = pos - 1
    if not (0 <= idx0 < len(sequence)) or sequence[idx0] != wt:
        raise ValueError(f"Position {pos}: expected {wt}, got {sequence[idx0] if idx0 < len(sequence) else 'OOB'}")

    tokens = tokenizer(
        sequence,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=1024,
    )
    input_ids = tokens["input_ids"].to(device)
    mask_idx = idx0 + 1  # +1 for <cls> token

    masked_ids = input_ids.clone()
    masked_ids[0, mask_idx] = tokenizer.mask_token_id

    with torch.no_grad():
        logits = model(masked_ids).logits

    log_probs = torch.log_softmax(logits[0, mask_idx], dim=-1)
    wt_id = tokenizer.convert_tokens_to_ids(wt)
    mt_id = tokenizer.convert_tokens_to_ids(mt)
    return (log_probs[mt_id] - log_probs[wt_id]).item()


def process_protein(fasta_path: str, var_path: str, model_names: list, device: torch.device) -> list:
    sequence = load_protein_sequence(fasta_path)
    variants = load_variants(var_path)
    if sequence is None or not variants:
        logger.warning("No sequence or variants in %s — skipping", fasta_path)
        return []

    results = []
    for model_name in model_names:
        model, tokenizer = load_model(model_name, device)
        for wt, pos, mt in variants:
            try:
                score = compute_variant_score(sequence, (wt, pos, mt), model, tokenizer, device)
                results.append({
                    "protein": os.path.basename(fasta_path),
                    "variant": f"{wt}{pos}{mt}",
                    "position": pos,
                    "ref_aa": wt,
                    "alt_aa": mt,
                    "score": score,
                    "model": model_name,
                })
            except Exception as exc:
                logger.error("Error scoring %s%d%s: %s", wt, pos, mt, exc)
    return results


# ---------------------------------------------------------------------------
# Analysis & output
# ---------------------------------------------------------------------------

def analyze_results(df: pd.DataFrame, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    plt.figure(figsize=(8, 4))
    sns.histplot(df["score"], bins=20, color="skyblue", edgecolor="black")
    plt.axvline(0, color="red", linestyle="--")
    plt.xlabel("Score (LLR)")
    plt.title("Distribution of Variant Scores")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "score_histogram.png"), dpi=150)
    plt.close()

    counts = [(df["score"] < 0).sum(), (df["score"] >= 0).sum()]
    plt.figure(figsize=(4, 4))
    plt.pie(counts, labels=["Deleterious", "Tolerated"], autopct="%1.1f%%",
            colors=["lightcoral", "lightgreen"])
    plt.title("Predicted Effect Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "effect_pie.png"), dpi=150)
    plt.close()

    if df["protein"].nunique() > 1:
        plt.figure(figsize=(8, 4))
        sns.boxplot(x="protein", y="score", data=df)
        plt.xticks(rotation=45)
        plt.title("Variant Scores by Protein")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "score_box_by_protein.png"), dpi=150)
        plt.close()

    logger.info("Figures written to %s", out_dir)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="ESM2 variant effect predictor")
    parser.add_argument("--fasta_dir", required=True, help="Directory with .fasta files")
    parser.add_argument("--var_dir", required=True, help="Directory with .var.txt files")
    parser.add_argument("--output", default="esm2_results.csv", help="Output CSV path")
    parser.add_argument(
        "--models",
        nargs="+",
        default=["facebook/esm2_t33_650M_UR50D"],
        help="One or more HuggingFace model identifiers",
    )
    parser.add_argument("--no_gpu", action="store_true", help="Force CPU execution")
    parser.add_argument("--figures_dir", default="figures", help="Directory for output plots")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    device = torch.device("cpu" if args.no_gpu or not torch.cuda.is_available() else "cuda")
    logger.info("Device: %s", device)
    if device.type == "cuda":
        logger.info("GPU: %s  (%.1f GB)", torch.cuda.get_device_name(), torch.cuda.get_device_properties(0).total_memory / 1e9)

    matches = find_matching_files(args.fasta_dir, args.var_dir)
    logger.info("Proteins to process: %d", len(matches))
    if not matches:
        logger.error("No matching fasta/variant pairs found — exiting")
        sys.exit(1)

    all_results = []
    for fasta_file, var_file, _ in tqdm(matches, desc="Proteins", unit="protein"):
        fasta_path = os.path.join(args.fasta_dir, fasta_file)
        all_results.extend(process_protein(fasta_path, var_file, args.models, device))

    if not all_results:
        logger.error("No results generated — exiting")
        sys.exit(1)

    df = pd.DataFrame(all_results)
    df["interpretation"] = df["score"].apply(lambda s: "Deleterious" if s < 0 else "Tolerated")

    logger.info("Proteins processed : %d", df["protein"].nunique())
    logger.info("Variants analysed  : %d", len(df))
    logger.info("Mean score         : %.3f", df["score"].mean())
    logger.info("Score range        : %.3f – %.3f", df["score"].min(), df["score"].max())

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, index=False)
    logger.info("Results saved to %s", args.output)

    analyze_results(df, args.figures_dir)


if __name__ == "__main__":
    main()
