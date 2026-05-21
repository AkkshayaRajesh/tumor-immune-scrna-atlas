"""08_composition.py — Cell type composition by EGFR status."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import mannwhitneyu
 
adata = sc.read_h5ad('data/processed/annotated.h5ad')
 
# --- Per-sample composition table ---
comp_per_sample = pd.crosstab(
    adata.obs['sample_id'], adata.obs['cell_type'],
    normalize='index') * 100
comp_per_sample.to_csv('results/tables/composition_per_sample.csv')
 
# --- Per-EGFR-group composition ---
comp_egfr = pd.crosstab(
    adata.obs['cell_type'], adata.obs['egfr_status'],
    normalize='columns') * 100
comp_egfr.to_csv('results/tables/composition_by_egfr.csv')
 
# --- Stacked bar plot ---
fig, ax = plt.subplots(figsize=(7, 6))
comp_egfr.T.plot.bar(stacked=True, colormap='tab20', ax=ax)
ax.set_ylabel('% of cells')
ax.set_xlabel('EGFR status')
ax.set_title('Cell type composition by EGFR status')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig('results/figures/composition_stacked.png', dpi=300)
plt.close()

# --- Statistical test per cell type, using patient-level proportions ---
# (per-cell tests inflate p-values - cells within a patient aren't independent)
results = []

for ct in adata.obs['cell_type'].unique():
    per_pt = comp_per_sample[ct]
    sample_egfr = adata.obs.drop_duplicates('sample_id').set_index(
        'sample_id')['egfr_status']
    mut = per_pt[sample_egfr.loc[per_pt.index] == 'mutant'].values
    wt  = per_pt[sample_egfr.loc[per_pt.index] == 'wildtype'].values
    if len(mut) >= 2 and len(wt) >= 2:
        stat, p = mannwhitneyu(mut, wt, alternative='two-sided')
        results.append({
            'cell_type':   ct,
            'mut_median':  np.median(mut),
            'wt_median':   np.median(wt),
            'log2_FC':     (np.median(mut) + 1) / (np.median(wt) + 1),
            'pvalue':      p,
        })
 
stats_df = pd.DataFrame(results).sort_values('pvalue')
stats_df.to_csv('results/tables/composition_stats.csv', index=False)
print(stats_df.to_string(index=False))
 
# --- Box plot of MDSC fraction by EGFR (the key panel) ---
if 'MDSC' in comp_per_sample.columns:
    plot_df = comp_per_sample[['MDSC']].copy()
    plot_df['egfr_status'] = sample_egfr.loc[plot_df.index].values
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.boxplot(data=plot_df, x='egfr_status', y='MDSC', ax=ax,
                hue='egfr_status', legend=False,
                palette={'mutant': '#d62728', 'wildtype': '#2ca02c'})
    sns.stripplot(data=plot_df, x='egfr_status', y='MDSC',
                  ax=ax, color='black', size=8)
    ax.set_title('MDSC fraction by EGFR status')
    ax.set_ylabel('% MDSC')
    plt.tight_layout()
    plt.savefig('results/figures/mdsc_by_egfr.png', dpi=300)
    plt.close()