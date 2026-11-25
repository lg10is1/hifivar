#!/usr/bin/env python3
"""
HiFiVar Benchmark Visualization Script

Generates publication-quality figures for HG002 benchmark results.
All data is derived from frozen, audited artifacts.
"""

import argparse
import json
from pathlib import Path

from benchmark_qualification import select_cmrg_callers

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Configuration
VISUAL_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = VISUAL_ROOT / 'figures' / 'benchmark'
DATA_DIR = VISUAL_ROOT / 'data'

MM_PER_INCH = 25.4
SINGLE_COLUMN = (89 / MM_PER_INCH, 72 / MM_PER_INCH)
DOUBLE_COLUMN = (183 / MM_PER_INCH, 96 / MM_PER_INCH)
DOUBLE_COLUMN_TALL = (183 / MM_PER_INCH, 116 / MM_PER_INCH)

# Colors
COLORS = {
    'precision': '#0072B2',
    'recall': '#009E73',
    'f1': '#E69F00',
    'gt': '#CC79A7',
    'tp': '#0072B2',
    'fp': '#D55E00',
    'fn': '#999999',
}

MARKERS = {'precision': 'o', 'recall': '^', 'f1': 's', 'gt': 'D'}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 6.5,
    'axes.titlesize': 7,
    'axes.labelsize': 6.5,
    'axes.linewidth': 0.6,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'legend.fontsize': 6,
    'legend.frameon': False,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'savefig.facecolor': 'white',
})


def style_axes(ax):
    """Apply the shared journal-style axis treatment."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(top=False, right=False)

def load_data(filename):
    """Load TSV data file."""
    filepath = DATA_DIR / filename
    return pd.read_csv(filepath, sep='\t')


def load_qualified_cmrg_data():
    """Exclude identity-unresolved results; never rewrite the archival TSV."""
    callers, _ = select_cmrg_callers(
        DATA_DIR / 'hg002_cmrg_sv_benchmark.tsv',
        DATA_DIR / 'cmrg_figure_qualifications.json',
    )
    df = load_data('hg002_cmrg_sv_benchmark.tsv')
    return df[df['caller'].isin(callers)].reset_index(drop=True)

def plot_small_variant_overview():
    """Plot 1: Small-variant performance overview."""
    df = load_data('hg002_small_variant_benchmark.tsv')
    
    fig, ax = plt.subplots(figsize=SINGLE_COLUMN)
    
    x = range(len(df))
    width = 0.25
    
    bars1 = ax.bar([i - width for i in x], df['precision'], width, 
                   label='Precision', color=COLORS['precision'])
    bars2 = ax.bar(x, df['recall'], width, 
                   label='Recall', color=COLORS['recall'])
    bars3 = ax.bar([i + width for i in x], df['f1'], width, 
                   label='F1', color=COLORS['f1'])
    
    ax.set_xlabel('Variant Class')
    ax.set_ylabel('Score')
    ax.set_title('HG002 small-variant performance', pad=24)
    ax.set_xticks(x)
    ax.set_xticklabels(df['variant_class'])
    ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    ax.set_ylim(0, 1.08)
    style_axes(ax)
    
    # Add value labels
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=5.5)
    
    plt.tight_layout()
    return fig

def plot_stratified_snp():
    """Plot 2: Stratified SNP performance."""
    df = load_data('hg002_small_variant_stratified.tsv')
    df_snp = df[df['variant_class'] == 'SNP'].copy()
    
    fig, ax = plt.subplots(figsize=DOUBLE_COLUMN_TALL)
    
    x = range(len(df_snp))
    width = 0.25
    
    bars1 = ax.bar([i - width for i in x], df_snp['precision'], width, 
                   label='Precision', color=COLORS['precision'])
    bars2 = ax.bar(x, df_snp['recall'], width, 
                   label='Recall', color=COLORS['recall'])
    bars3 = ax.bar([i + width for i in x], df_snp['f1'], width, 
                   label='F1', color=COLORS['f1'])
    
    ax.set_xlabel('Stratification')
    ax.set_ylabel('Score')
    ax.set_title('HG002 stratified SNP performance')
    ax.set_xticks(x)
    labels = df_snp['stratification'].str.replace('_', ' ', regex=False)
    ax.set_xticklabels(labels, rotation=35, ha='right')
    ax.legend(ncol=3)
    ax.set_ylim(0, 1.05)
    style_axes(ax)
    
    plt.tight_layout()
    return fig

def plot_stratified_indel():
    """Plot 3: Stratified INDEL performance."""
    df = load_data('hg002_small_variant_stratified.tsv')
    df_indel = df[df['variant_class'] == 'INDEL'].copy()
    
    fig, ax = plt.subplots(figsize=DOUBLE_COLUMN_TALL)
    
    x = range(len(df_indel))
    width = 0.25
    
    bars1 = ax.bar([i - width for i in x], df_indel['precision'], width, 
                   label='Precision', color=COLORS['precision'])
    bars2 = ax.bar(x, df_indel['recall'], width, 
                   label='Recall', color=COLORS['recall'])
    bars3 = ax.bar([i + width for i in x], df_indel['f1'], width, 
                   label='F1', color=COLORS['f1'])
    
    ax.set_xlabel('Stratification')
    ax.set_ylabel('Score')
    ax.set_title('HG002 stratified INDEL performance')
    ax.set_xticks(x)
    labels = df_indel['stratification'].str.replace('_', ' ', regex=False)
    ax.set_xticklabels(labels, rotation=35, ha='right')
    ax.legend(ncol=3)
    ax.set_ylim(0, 1.05)
    style_axes(ax)
    
    plt.tight_layout()
    return fig

def plot_cmrg_sv_comparison():
    """Plot 4: CMRG SV caller comparison."""
    df = load_qualified_cmrg_data()
    
    fig, ax = plt.subplots(figsize=DOUBLE_COLUMN)
    
    x = range(len(df))
    width = 0.2
    
    bars1 = ax.bar([i - 1.5*width for i in x], df['precision'], width, 
                   label='Precision', color=COLORS['precision'])
    bars2 = ax.bar([i - 0.5*width for i in x], df['recall'], width, 
                   label='Recall', color=COLORS['recall'])
    bars3 = ax.bar([i + 0.5*width for i in x], df['f1'], width, 
                   label='F1', color=COLORS['f1'])
    bars4 = ax.bar([i + 1.5*width for i in x], df['gt_concordance'], width, 
                   label='GT concordance', color=COLORS['gt'])
    
    ax.set_xlabel('SV Caller')
    ax.set_ylabel('Score')
    ax.set_title('HG002 CMRG historical SV estimates', pad=24)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{row['caller']}\n{row['version']}" for _, row in df.iterrows()])
    ax.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    ax.set_ylim(0.0, 1.05)
    style_axes(ax)
    
    plt.tight_layout()
    return fig

def plot_cmrg_sv_breakdown():
    """Plot 5: CMRG SV TP/FP/FN breakdown."""
    df = load_qualified_cmrg_data()
    
    fig, ax = plt.subplots(figsize=DOUBLE_COLUMN)
    
    x = range(len(df))
    width = 0.6
    
    bars1 = ax.bar(x, df['tp'], width, label='TP', color=COLORS['tp'])
    bars2 = ax.bar(x, df['fp'], width, bottom=df['tp'], label='FP', color=COLORS['fp'])
    bars3 = ax.bar(x, df['fn'], width, bottom=df['tp'] + df['fp'], label='FN', color=COLORS['fn'])
    
    ax.set_xlabel('SV Caller')
    ax.set_ylabel('Count')
    ax.set_title('HG002 CMRG historical SV outcomes', pad=24)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{row['caller']}\n{row['version']}" for _, row in df.iterrows()])
    ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    style_axes(ax)
    
    plt.tight_layout()
    return fig

def plot_precision_recall():
    """Plot 6: Small-variant precision-recall."""
    df = load_data('hg002_small_variant_benchmark.tsv')
    
    fig, ax = plt.subplots(figsize=SINGLE_COLUMN)
    
    colors = [COLORS['precision'], COLORS['f1']]
    markers = ['o', 's']
    
    for i, row in df.iterrows():
        ax.scatter(row['recall'], row['precision'], 
                   s=36, c=colors[i], marker=markers[i],
                   label=row['variant_class'], edgecolors='black', linewidths=0.5)
    
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('HG002 small-variant precision–recall')
    ax.legend()
    ax.set_xlim(0.85, 1.0)
    ax.set_ylim(0.85, 1.0)
    ax.set_aspect('equal')
    style_axes(ax)
    
    plt.tight_layout()
    return fig

def save_figure(fig, filename):
    """Save figure in multiple formats."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    base = OUTPUT_DIR / filename
    fig.savefig(f'{base}.svg', format='svg', bbox_inches='tight')
    fig.savefig(f'{base}.pdf', format='pdf', bbox_inches='tight')
    fig.savefig(f'{base}.png', format='png', dpi=300, bbox_inches='tight')
    
    plt.close(fig)
    print(f"Saved: {filename}.svg, {filename}.pdf, {filename}.png")

def main():
    """Generate all benchmark figures."""
    global OUTPUT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    OUTPUT_DIR = args.output_dir.resolve()
    _, qualification = select_cmrg_callers(
        DATA_DIR / 'hg002_cmrg_sv_benchmark.tsv',
        DATA_DIR / 'cmrg_figure_qualifications.json',
    )
    print("Generating HiFiVar benchmark figures...")
    
    # Figure 1
    fig = plot_small_variant_overview()
    save_figure(fig, 'hifivar_benchmark_small_variant_overview')
    
    # Figure 2
    fig = plot_stratified_snp()
    save_figure(fig, 'hifivar_benchmark_stratified_snp')
    
    # Figure 3
    fig = plot_stratified_indel()
    save_figure(fig, 'hifivar_benchmark_stratified_indel')
    
    # Figure 4
    fig = plot_cmrg_sv_comparison()
    save_figure(fig, 'hifivar_benchmark_cmrg_sv_comparison')
    
    # Figure 5
    fig = plot_cmrg_sv_breakdown()
    save_figure(fig, 'hifivar_benchmark_cmrg_sv_breakdown')
    
    # Figure 6
    fig = plot_precision_recall()
    save_figure(fig, 'hifivar_benchmark_precision_recall')
    
    # Revalidate source identity after plotting; no changed input can be frozen.
    _, final_qualification = select_cmrg_callers(
        DATA_DIR / 'hg002_cmrg_sv_benchmark.tsv',
        DATA_DIR / 'cmrg_figure_qualifications.json',
    )
    if final_qualification != qualification:
        raise ValueError('CMRG figure qualification changed during generation')
    (OUTPUT_DIR / 'figure_qualification.json').write_text(
        json.dumps(qualification, indent=2) + '\n', encoding='utf-8',
    )
    print("All figures generated successfully; qualification manifest retained.")

if __name__ == '__main__':
    main()
