"""05a_harmony.py — Harmony integration."""
import os
import scanpy as sc
import harmonypy as hm
 
IN_PATH  = 'data/interim/04_normalized.h5ad'
OUT_PATH = 'data/interim/05a_harmony.h5ad'
 
adata = sc.read_h5ad(IN_PATH)
 
# Scale (Harmony works on PCs, which work better on scaled data)
sc.pp.scale(adata,max_value=10)

#PCA
sc.tl.pca(adata,n_comps=50)

#Harmony — call directly to avoid scanpy wrapper double-transposition bug
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample_id')
adata.obsm['X_pca_harmony'] = ho.Z_corr  # this harmonypy version returns Z_corr as (n_cells, n_dims)

# Neighbors + UMAP on Harmony embedding
sc.pp.neighbors(adata,use_rep='X_pca_harmony', n_neighbors=15)
sc.tl.umap(adata)
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
adata.write(OUT_PATH)
print(f'Wrote {OUT_PATH}')