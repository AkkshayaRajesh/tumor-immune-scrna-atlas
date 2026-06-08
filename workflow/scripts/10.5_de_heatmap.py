import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os

adata = sc.read_h5ad('data/processed/annotated.h5ad')

top_genes = {}
for ct in adata.obs['cell_type'].unique():
    f = f'results/tables/de_pseudobulk/de_{ct}.csv'
    if not os.path.exists(f): continue
    df = pd.read_csv(f)
    df = df[df['padj'] < 0.05].sort_values(
        'log2FoldChange', ascending=False)
    top_genes[ct] = df.head(5)['gene'].tolist()
 
# Heatmap: scaled expression of top genes across cell types
all_genes = sorted({g for genes in top_genes.values() for g in genes})
if not all_genes:
    raise ValueError('No significant DE genes found — check de_pseudobulk CSVs exist and have padj < 0.05 hits')

sc.pl.heatmap(adata, var_names=all_genes,
              groupby=['cell_type', 'egfr_status'],
              standard_scale='var', dendrogram=False,
              save='_de_top_genes.png')