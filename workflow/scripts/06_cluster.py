"""06_cluster.py — Multi-resolution Leiden clustering."""
import os
import scanpy as sc
 
IN_PATH  = 'data/interim/05_integrated.h5ad'
OUT_PATH = 'data/interim/06_clustered.h5ad'
 
adata = sc.read_h5ad(IN_PATH)

# Make sure we're using scVI neighbors graph
sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=15)

for res in [0.3, 0.6, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    n = adata.obs[f'leiden_{res}'].nunique()
    print(f'Resolution {res}: {n} clusters')

# Panel figure of all four
sc.pl.umap(
    adata, color=['leiden_0.3', 'leiden_0.6',
                  'leiden_1.0', 'leiden_1.5'],
    ncols=2, save='_leiden_resolutions.png', legend_loc='on data',
    legend_fontsize=6,
)
 
# Set the working resolution — adjust after inspecting the figure
WORKING_RES = '0.6'   # changed after visual inspection
''' "I chose resolution 0.6 (18 clusters) because it gives clear T cell, NK, B, plasma, myeloid, and epithelial separation without over-fragmenting the myeloid compartment."'''

adata.obs['leiden'] = adata.obs[f'leiden_{WORKING_RES}']
print(f'Working resolution: {WORKING_RES}, ',
      f'{adata.obs["leiden"].nunique()} clusters')
 
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
adata.write(OUT_PATH)