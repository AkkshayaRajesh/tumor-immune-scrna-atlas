"""02_qc.py — QC filtering and metric computation."""

import os
import scanpy as sc

try:
    IN_PATH  = snakemake.input[0]    # noqa
    OUT_PATH = snakemake.output[0]   # noqa
except NameError:
    IN_PATH  = 'data/raw/01_concatenated.h5ad'
    OUT_PATH = 'data/interim/02_qc.h5ad'

sc.settings.figdir = 'results/figures/'
sc.set_figure_params(dpi_save=900)


adata = sc.read_h5ad(IN_PATH)

# Mitochondrial and ribosomal annotations
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS','RPL'))

sc.pp.calculate_qc_metrics(
    adata,qc_vars=['mt','ribo'],
    inplace=True,percent_top=None,log1p=False
)

#BEFORE plots
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='sample_id', rotation=45,multi_panel=True,
    save='_qc_before.png',
)

#Filtering
n_before = adata.n_obs
sc.pp.filter_cells(adata,min_genes=200)
sc.pp.filter_genes(adata,min_cells=3)
adata = adata[adata.obs['n_genes_by_counts']<6000].copy()
adata = adata[adata.obs['pct_counts_mt']<20].copy()
n_after = adata.n_obs

print(f'Filtered: {n_before:,} -> {n_after:,}'
      f'(removed {n_before-n_after:,},'
      f'{100*(n_before-n_after)/n_before:.1f}%)')

#AFTER plots
sc.pl.violin(
    adata,
    ['n_genes_by_counts','total_counts','pct_counts_mt'],
    groupby='sample_id', rotation=45, multi_panel=True,
    save='_qc_after.png'
)


#Per-sample yield table
yields = adata.obs.groupby('sample_id').size().reset_index(
    name='cells_post_qc')
yields.to_csv('results/tables/qc_yields.csv', index=False)
 
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
adata.write(OUT_PATH)
print(f'Wrote {OUT_PATH}')