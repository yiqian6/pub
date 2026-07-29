#!/usr/bin/env python3
"""
Overlap TRACE windows with gene annotations (UCSC refGene) and filter for target genes.

Usage:
    python gene_filter.py --report_dir /path/to/trace/report/ [--target_genes OSM,LIF] [--output overlap.csv] [--gene_annot refGene.txt]

Arguments:
    --report_dir    : TRACE report directory containing trace_all_gw.csv (required)
    --target_genes  : Comma-separated gene symbols (default: OSM,LIF,GAL3ST1,NF2)
    --output        : Output CSV file (default: trace_gene_overlap.csv)
    --gene_annot    : Local refGene.txt file (if not provided, download from UCSC)
"""

import argparse
import os
import gzip
import urllib.request
import shutil
import pandas as pd

def get_gene_bed(gene_file='refGene.txt'):
    """Download and parse UCSC refGene into BED-like DataFrame (chr22 only)."""
    if not os.path.exists(gene_file):
        print("Downloading refGene.txt.gz from UCSC...")
        url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"
        try:
            urllib.request.urlretrieve(url, 'refGene.txt.gz')
            with gzip.open('refGene.txt.gz', 'rb') as f_in:
                with open(gene_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove('refGene.txt.gz')
        except Exception as e:
            raise SystemExit(f"Failed to download gene annotation: {e}")

    print("Parsing gene annotation...")
    df_all = pd.read_csv(gene_file, sep='\t', header=None)
    # Column indices: 2=chrom, 4=txStart, 5=txEnd, 12=name2
    df = df_all[[2,4,5,12]].copy()
    df.columns = ['chrom', 'start', 'end', 'gene']
    df['chrom'] = df['chrom'].str.replace('chr', '').astype(str)
    df = df[df['chrom'] == '22']
    df = df.groupby('gene').agg({
        'chrom': 'first',
        'start': 'min',
        'end': 'max'
    }).reset_index()
    df[['start','end']] = df[['start','end']].astype(int)
    print(f"Loaded {len(df)} genes on chr22")
    return df

def get_trace_windows(report_dir):
    gw_file = os.path.join(report_dir, 'trace_all_gw.csv')
    if not os.path.isfile(gw_file):
        raise SystemExit(f"Error: trace_all_gw.csv not found in {report_dir}")
    df = pd.read_csv(gw_file)
    df['chrom'] = df['chrom'].astype(str)
    df['start_0'] = df['win_start'] - 1
    df['end_0'] = df['win_end']
    print(f"Loaded {len(df)} windows from {gw_file}")
    return df

def find_all_overlaps(trace_df, gene_df):
    overlaps = []
    for _, t_row in trace_df.iterrows():
        chrom = t_row['chrom']
        s0 = t_row['start_0']
        e0 = t_row['end_0']
        matched = gene_df[
            (gene_df['chrom'] == chrom) &
            (gene_df['start'] < e0) &
            (gene_df['end'] > s0)
        ]
        if not matched.empty:
            gene_names = matched['gene'].unique()
            overlaps.append({
                'window': t_row['window'],
                'win_start': t_row['win_start'],
                'win_end': t_row['win_end'],
                'burden_percent': t_row['burden_percent'],
                'n_samples': t_row['n_samples_with_segment'],
                'genes': ','.join(gene_names),
                'gene_count': len(gene_names)
            })
    return pd.DataFrame(overlaps)

def main():
    parser = argparse.ArgumentParser(description="Filter TRACE windows by target genes")
    parser.add_argument("--report_dir", required=True, help="TRACE report directory")
    parser.add_argument("--target_genes", default="OSM,LIF,GAL3ST1,NF2",
                        help="Comma-separated gene symbols (default: OSM,LIF,GAL3ST1,NF2)")
    parser.add_argument("--output", default="trace_gene_overlap.csv",
                        help="Output CSV file (default: trace_gene_overlap.csv)")
    parser.add_argument("--gene_annot", default="refGene.txt",
                        help="Local refGene.txt file (if not found, download)")
    args = parser.parse_args()

    target_genes = [g.strip() for g in args.target_genes.split(',') if g.strip()]

    # Get gene annotation (download if needed)
    gene_df = get_gene_bed(args.gene_annot)

    # Filter to target genes
    target_gene_df = gene_df[gene_df['gene'].isin(target_genes)]
    if target_gene_df.empty:
        print(f"Warning: None of the target genes {target_genes} found in annotation.")
        print("Available genes on chr22 (sample):", gene_df['gene'].head(10).tolist())
        return

    # Load TRACE windows
    trace_df = get_trace_windows(args.report_dir)

    # Find overlaps
    result = find_all_overlaps(trace_df, target_gene_df)
    if result.empty:
        print("No windows overlap with any target genes.")
        return

    result = result.sort_values('burden_percent', ascending=False)
    result.to_csv(args.output, index=False)
    print(f"Found {len(result)} windows covering target genes.")
    print(f"Results saved to {args.output}")
    print("\nTop 10 windows by burden:")
    print(result[['window', 'burden_percent', 'genes']].head(10))

if __name__ == "__main__":
    main()