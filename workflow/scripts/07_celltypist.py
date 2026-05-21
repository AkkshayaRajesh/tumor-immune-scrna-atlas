"""07_celltypist.py — Automated annotation with CellTypist."""
import scanpy as sc
import celltypist
from celltypist import models


IN_PATH  = 'data/interim/06_clustered.h5ad'
OUT_PATH = 'data/interim/07_celltypist.h5ad'

adata = sc.read_h5ad(IN_PATH)

# Download the immune cell model — well-suited for TIL data
models.download_models(model='Immune_All_Low.pkl')

# CellTypist expects log-normalized counts with 10k normalization
# adata.X is already that (HVG-subset), but for annotation we want he full gene set in adata.raw
adata_full = adata.raw.to_adata()
adata_full.obs = adata.obs.copy()

predictions= celltypist.annotate(
    adata_full,
    model='Immune_All_Low.pkl',
    majority_voting=True,
    over_clustering='leiden',
)

#Transfer predictions to our HVG-subset object
adata.obs['celltypist_label'] = predictions.predicted_labels['majority_voting'].astype('category')
adata.obs['celltypist_over_clustering'] = predictions.predicted_labels[
    'over_clustering'].astype('category')


sc.pl.umap(adata,color='celltypist_label',
           legend_loc='on data', legend_fontsize=6,
           save='_celltypist.png')

adata.write(OUT_PATH)