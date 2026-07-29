#!/usr/bin/env python3
"""
Generate burden line plot and carrier histogram from TRACE output.

Usage:
    python visualize_trace.py --report_dir /path/to/trace/report/ [--output_dir ./figures]

Arguments:
    --report_dir   : TRACE report directory containing trace_all_gw.csv (required)
    --output_dir   : Directory to save figures (default: ./figures)
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Visualize TRACE burden results")
    parser.add_argument("--report_dir", required=True, help="TRACE report directory")
    parser.add_argument("--output_dir", default="./figures", help="Output directory for figures")
    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # Read data
    gw_file = os.path.join(args.report_dir, "trace_all_gw.csv")
    if not os.path.isfile(gw_file):
        raise SystemExit(f"Error: trace_all_gw.csv not found in {args.report_dir}")

    df = pd.read_csv(gw_file)
    df['pos_mb'] = df['win_start'] / 1_000_000
    chr_name = df['chrom'].iloc[0]

    print(f"Total windows: {len(df)}, Chromosome: {chr_name}")
    print(f"Max burden: {df['burden_percent'].max():.2f}% at window {df.loc[df['burden_percent'].idxmax(), 'window']}")

    # ----------------------------------------
    # Figure 1: Burden line plot
    # ----------------------------------------
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.plot(df['pos_mb'], df['burden_percent'], color='#2E86C1', linewidth=1.8, label='Burden (%)')
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.4, label='100% baseline')

    max_row = df.loc[df['burden_percent'].idxmax()]
    ax.annotate(f'Peak {max_row["burden_percent"]:.1f}%',
                xy=(max_row['pos_mb'], max_row['burden_percent']),
                xytext=(max_row['pos_mb'] + 1.5, max_row['burden_percent'] - 30),
                arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5),
                fontsize=11, color='darkred')

    ax.set_xlabel(f'Position on {chr_name} (Mb)', fontsize=13)
    ax.set_ylabel('Cumulative sample burden (%)', fontsize=13)
    ax.set_title(f'TRACE_Archaic segment prevalence on {chr_name}', fontsize=15)
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right')
    plt.tight_layout()
    line_plot = os.path.join(args.output_dir, 'burden_line.png')
    plt.savefig(line_plot, dpi=300)
    plt.close()
    print(f"Burden line plot saved to {line_plot}")

    # ----------------------------------------
    # Figure 2: Carrier histogram
    # ----------------------------------------
    fig2, ax2 = plt.subplots(figsize=(10, 4))
    ax2.hist(df['n_samples_with_segment'], bins=25, color='#28B463', edgecolor='black', alpha=0.7)
    ax2.axvline(x=df['n_samples_with_segment'].median(),
                color='blue', linestyle='--',
                label=f"Median = {df['n_samples_with_segment'].median():.0f}")
    ax2.set_xlabel('Number of samples carrying segment per 1Mb window', fontsize=12)
    ax2.set_ylabel('Number of windows', fontsize=12)
    ax2.set_title('Distribution of segment carriers across chr22', fontsize=14)
    ax2.legend()
    plt.tight_layout()
    hist_plot = os.path.join(args.output_dir, 'carrier_hist.png')
    plt.savefig(hist_plot, dpi=300)
    plt.close()
    print(f"Carrier histogram saved to {hist_plot}")

if __name__ == "__main__":
    main()