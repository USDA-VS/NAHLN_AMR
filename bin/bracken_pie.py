#!/usr/bin/env python3
"""
Generate Bracken pie chart from Excel file
Usage: bracken_pie.py sample_name bracken.xlsx
Written by RW and Gemini 2.5 Pro under USAi 2026
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')

def generate_pie_chart(sample_name, bracken_excel, output_png=None):
    """Generate pie chart from Bracken Excel file"""
    
    try:
        df = pd.read_excel(bracken_excel)
        if 'fraction_total_reads' not in df.columns or df.empty:
            raise ValueError("Input Excel file is empty or missing 'fraction_total_reads' column.")
    except Exception as e:
        print(f"ERROR: Could not read or parse {bracken_excel}. {e}")
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.text(0.5, 0.5, f'Error processing Bracken results\nfor {sample_name}',
                ha='center', va='center', fontsize=14, color='red',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')
        if output_png is None:
            output_png = f"{sample_name}_bracken_pie.png"
        fig.savefig(output_png, format='png', bbox_inches='tight', dpi=150)
        return output_png
    
    numeric_cols = df.select_dtypes(include=['number']).columns
    df[numeric_cols] = df[numeric_cols].clip(lower=0)
    df = df[df['fraction_total_reads'] > 0.01]

    if len(df) == 0:
        print(f"WARNING: No species with >1% abundance for {sample_name}")
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.text(0.5, 0.5, f'Insufficient classified reads\nfor species-level identification\n({sample_name})',
                ha='center', va='center', fontsize=14,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')
        if output_png is None: output_png = f"{sample_name}_bracken_pie.png"
        fig.savefig(output_png, format='png', bbox_inches='tight', dpi=150)
        return output_png

    classified_fraction = df['fraction_total_reads'].sum()
    unclassified_fraction = max(0, 1 - classified_fraction)

    threshold = 0.02
    df_main = df[df['fraction_total_reads'] >= threshold]
    other_fraction = df[df['fraction_total_reads'] < threshold]['fraction_total_reads'].sum()

    if other_fraction > 0.001:
        df_other = pd.DataFrame([{'name': f'Other (<{threshold:.0%})', 'fraction_total_reads': other_fraction}])
        df_processed = pd.concat([df_main, df_other])
    else:
        df_processed = df_main

    if unclassified_fraction > 0.001:
        df_unclassified = pd.DataFrame([{'name': 'unclassified', 'fraction_total_reads': unclassified_fraction}])
        df_processed = pd.concat([df_processed, df_unclassified])

    df_processed = df_processed.set_index('name')
    
    cmap = plt.get_cmap('coolwarm')
    fig, ax = plt.subplots(figsize=(9, 5), dpi=150)
    
    df_processed.plot.pie(
        y='fraction_total_reads',
        ax=ax,
        title='FASTQ Read Identification',
        cmap=cmap,
        labeldistance=None,
        legend=True,
        autopct=lambda p: f'{p:.1f}%' if p > 2 else '',
        pctdistance=0.7,
        textprops={'color': 'black', 'fontsize': 10}
    )
    
    ax.axis('off')
    ax.legend(bbox_to_anchor=(0.9, 0.9))
    ax.yaxis.label.set_visible(False)
    
    if other_fraction > 0.001:
        note = f"*Note: Species with < {threshold:.0%} abundance are grouped into the 'Other' category."
        fig.text(0.5, 0.05, note, ha='center', va='center', fontsize=8, fontstyle='italic', color='gray')
    
    if output_png is None:
        output_png = f"{sample_name}_bracken_pie.png"
    fig.savefig(output_png, format='png', bbox_inches='tight')
    print(f"Pie chart saved: {output_png}")
    return output_png

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: bracken_pie.py sample_name bracken.xlsx [output.png]")
        sys.exit(1)
    sample_name = sys.argv[1]
    bracken_excel = sys.argv[2]
    output_png = sys.argv[3] if len(sys.argv) > 3 else None
    generate_pie_chart(sample_name, bracken_excel, output_png)