"""05b_scvi.py — scVI integration trained from scratch."""
import os
import scanpy as sc
import scvi

IN_PATH  = 'data/interim/05a_harmony.h5ad'
OUT_PATH = 'data/interim/05b_integrated.h5ad'

adata = sc.read_h5ad(IN_PATH)

# Raw counts stored in layers['counts'] by step 04 — required by scVI
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='sample_id')

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
model.train(max_epochs=400, early_stopping=True)

adata.obsm['X_scVI'] = model.get_latent_representation()

# Neighbors + UMAP on scVI embedding
sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=15, key_added='scvi')
sc.tl.umap(adata, neighbors_key='scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
adata.write(OUT_PATH)
print(f'Wrote {OUT_PATH}')

if os.path.exists('data/interim/05a_harmony.h5ad'):
    os.remove('data/interim/05a_harmony.h5ad')
