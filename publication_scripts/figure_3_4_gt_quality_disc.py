#!/usr/bin/env python3
"""
Analyse genotype concordance from a GATK GCsS file and plot:
  1. discordance_vs_quality.png  — scatter: depth / GQ vs discordance rate
  2. concordance_two_panel.png   — ranked discordance overview + labelled outliers

Usage:
    figure_3_gt_disc.py <gcss_file> [vcf_file] [output_dir]
"""

import os
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import vcf

def parse_gcss_file(filepath):
    """Parse GCsS and flag outliers using 2 standard deviations method"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    data_rows = []
    for line in lines:
        if line.startswith('GCsS') and not line.startswith('GCsS\t['):
            parts = line.strip().split('\t')
            if len(parts) >= 11:
                non_reference_discordance = float(parts[3].replace(',', '.'))

                # Calculate discordance for all three genotype classes
                rr_total = int(parts[4]) + int(parts[7])
                ra_total = int(parts[5]) + int(parts[8])
                aa_total = int(parts[6]) + int(parts[9])
                
                rr_concordance = 100 * int(parts[4]) / rr_total if rr_total > 0 else 0
                ra_concordance = 100 * int(parts[5]) / ra_total if ra_total > 0 else 0
                aa_concordance = 100 * int(parts[6]) / aa_total if aa_total > 0 else 0
                
                data_rows.append({
                    'sample': parts[2],
                    'discordance_rate': non_reference_discordance,
                    'rr_discordance': 100 - rr_concordance,
                    'ra_discordance': 100 - ra_concordance,
                    'aa_discordance': 100 - aa_concordance,
                })
    
    df = pd.DataFrame(data_rows)
    mean_discord = df['discordance_rate'].mean()
    std_discord = df['discordance_rate'].std()
    threshold = mean_discord + 2 * std_discord
    df['is_outlier'] = df['discordance_rate'] > threshold
    df['threshold'] = threshold  # Store threshold in dataframe
    df['discordance_mean'] = mean_discord
    df['discordance_std'] = std_discord
    
    return df

def extract_quality_metrics(vcf_file, gcss_df):
    """Extract quality metrics from VCF"""
    
    reader = vcf.Reader(filename=vcf_file)
    sample_names = set(gcss_df['sample'].values)
    
    metrics = defaultdict(lambda: {
        'dp_values': [],
        'gq_values': [],
        'ad_ref': [],
        'ad_alt': [],
    })
    
    print("Extracting quality metrics from VCF...")
    for i, record in enumerate(reader, 1):
        if i % 10000 == 0:
            print(f"  Processed {i} variants...")
        
        for sample in record.samples:
            if sample.sample not in sample_names:
                continue
            
            if sample['DP'] is not None:
                metrics[sample.sample]['dp_values'].append(sample['DP'])
            if sample['GQ'] is not None:
                metrics[sample.sample]['gq_values'].append(sample['GQ'])
            if sample['AD'] is not None and len(sample['AD']) >= 2:
                metrics[sample.sample]['ad_ref'].append(sample['AD'][0])
                metrics[sample.sample]['ad_alt'].append(sample['AD'][1])
    
    summary = []
    for sample_name, data in metrics.items():
        summary.append({
            'sample': sample_name,
            'mean_dp': np.mean(data['dp_values']) if data['dp_values'] else 0,
            'median_dp': np.median(data['dp_values']) if data['dp_values'] else 0,
            'std_dp': np.std(data['dp_values']) if data['dp_values'] else 0,
            'mean_gq': np.mean(data['gq_values']) if data['gq_values'] else 0,
            'median_gq': np.median(data['gq_values']) if data['gq_values'] else 0,
            'std_gq': np.std(data['gq_values']) if data['gq_values'] else 0,
        })
    
    return pd.DataFrame(summary)

def plot_figures(gcss_df, metrics_df, output_dir='./'):
    """Save discordance_vs_quality and concordance_two_panel figures."""
    os.makedirs(output_dir, exist_ok=True)

    df        = gcss_df.merge(metrics_df, on='sample', how='left')
    outliers  = df[df['is_outlier']]
    normal    = df[~df['is_outlier']]
    threshold = df['threshold'].iloc[0]

    # discordance_vs_quality
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Depth vs Discordance
    axes[0].scatter(normal['mean_dp'], normal['discordance_rate'], 
                   s=100, alpha=0.6, color='#0173B2', label='Normal', edgecolors='black')
    axes[0].scatter(outliers['mean_dp'], outliers['discordance_rate'], 
                   s=150, alpha=0.8, color='#DE8F05', label='High-Discordance', 
                   edgecolors='black', linewidth=2)
    
    axes[0].set_xlabel('Mean Depth (DP)', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Non-Reference Discordance (%)', fontsize=11, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)
    
    # GQ vs Discordance
    axes[1].scatter(normal['mean_gq'], normal['discordance_rate'], 
                   s=100, alpha=0.6, color='#0173B2', label='Normal', edgecolors='black')
    axes[1].scatter(outliers['mean_gq'], outliers['discordance_rate'], 
                   s=150, alpha=0.8, color='#DE8F05', label='High-Discordance', 
                   edgecolors='black', linewidth=2)
    
    axes[1].set_xlabel('Mean Genotype Quality (GQ)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Non-Reference Discordance (%)', fontsize=11, fontweight='bold')
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)
    
    plt.suptitle('Relationship Between Quality Metrics and Discordance', 
                fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/discordance_vs_quality.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: discordance_vs_quality.png")

    # concordance_two_panel
    mean_all = df['discordance_rate'].mean()
    std_all = df['discordance_rate'].std()
    mean_normal = df[~df['is_outlier']]['discordance_rate'].mean()
    std_normal = df[~df['is_outlier']]['discordance_rate'].std()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6),
                                   gridspec_kw={"width_ratios": [3, 1]})
    
    # Sort all samples by discordance for left panel
    df_sorted = df.sort_values('discordance_rate', ascending=True).reset_index(drop=True)
    
    # --- Left panel: Full overview ---
    colors = ['#DE8F05' if v >= threshold else '#0173B2'
              for v in df_sorted['discordance_rate']]
    ax1.bar(range(len(df_sorted)), df_sorted['discordance_rate'], 
            color=colors, width=1.0, linewidth=0, alpha=0.85)
    ax1.axhline(threshold, color='red', linestyle='--', linewidth=1.5, 
                label=f'Threshold ({threshold:.1f}%)')
    ax1.set_xticks([0, len(df_sorted) - 1])
    ax1.set_xticklabels(['Rank 1\n(Lowest)', f'Rank {len(df_sorted)}\n(Highest)'])
    ax1.set_ylabel('Non-Reference Discordance (%)', fontsize=12, fontweight='bold')
    ax1.set_xlabel(f'All samples (n={len(df_sorted)})', fontsize=12, fontweight='bold')
    ax1.set_title('Overview', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper left')
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_axisbelow(True)
    
    # Add mean concordance boxes on the left panel
    ax1.text(0.02, 0.50, f'Mean Discordance (All):\n{mean_all:.2f}%\nSD: {std_all:.2f}%', 
            transform=ax1.transAxes, fontsize=10, verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax1.text(0.02, 0.35, f'Mean Discordance (Normal):\n{mean_normal:.2f}%\nSD: {std_normal:.2f}%', 
            transform=ax1.transAxes, fontsize=10, verticalalignment='top',
            horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    # --- Right panel: Outliers only with sample labels ---
    outliers_sorted = df[df['is_outlier']].sort_values('discordance_rate', ascending=True)
    ax2.barh(range(len(outliers_sorted)), outliers_sorted['discordance_rate'],
             color='#DE8F05', linewidth=0.5, edgecolor='black', alpha=0.85)
    ax2.axvline(threshold, color='red', linestyle='--', linewidth=1.5,
                label=f'Threshold ({threshold:.1f}%)')
    ax2.set_yticks(range(len(outliers_sorted)))
    ax2.set_yticklabels(outliers_sorted['sample'], fontsize=9)
    ax2.set_xlabel('Non-Reference Discordance (%)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Outliers (n={len(outliers_sorted)})', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10, loc='lower right')
    ax2.spines[['top', 'right']].set_visible(False)
    ax2.grid(True, alpha=0.3, axis='x')
    ax2.set_axisbelow(True)
    
    plt.suptitle('Per-Sample Genotype Concordance', fontsize=15, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/concordance_two_panel.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: concordance_two_panel.png")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot genotype discordance figures from a GATK GCsS file.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  figure_3_gt_disc.py gcss_file.txt variants.vcf.gz ./output/'
    )
    parser.add_argument('gcss_file',
                        help='GCsS (Genotype Concordance Summary Statistics) file from GATK')
    parser.add_argument('vcf_file', nargs='?', default='0003.vcf.gz',
                        help='VCF file for quality metrics (default: 0003.vcf.gz)')
    parser.add_argument('output_dir', nargs='?', default='./',
                        help='Output directory (default: current directory)')
    args = parser.parse_args()

    gcss_df    = parse_gcss_file(args.gcss_file)
    metrics_df = extract_quality_metrics(args.vcf_file, gcss_df)
    plot_figures(gcss_df, metrics_df, args.output_dir)