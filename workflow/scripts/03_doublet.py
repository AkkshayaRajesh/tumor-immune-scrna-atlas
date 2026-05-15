
"""03_doublet.py — Per-sample doublet detection with Scrublet."""
import os
import numpy as np
import scanpy as sc

try:
    IN_PATH  = snakemake.input[0]    # noqa
    OUT_PATH = snakemake.output[0]   # noqa
except NameError:
    IN_PATH  = 'data/interim/02_qc.h5ad'
    OUT_PATH = 'data/interim/03_doublet.h5ad'

adata = sc.read_h5ad(IN_PATH)
adata.obs['doublet_score']     = np.nan
adata.obs['predicted_doublet'] = False

for sid in adata.obs['sample_id'].unique():
    mask = adata.obs['sample_id'] == sid
    sub = adata[mask].copy()
    if sub.n_obs < 100:
        print(f'{sid}: skipping ({sub.n_obs} cells)')
        continue
    try:
        sc.external.pp.scrublet(sub, verbose=False)
        scores    = sub.obs['doublet_score'].values
        predicted = sub.obs['predicted_doublet'].values
        adata.obs.loc[mask, 'doublet_score']     = scores
        adata.obs.loc[mask, 'predicted_doublet'] = predicted
        print(f'{sid}: {predicted.sum():,} / {sub.n_obs:,} doublets '
              f'({100*predicted.mean():.1f}%)')
    except Exception as e:
        print(f'{sid}: Scrublet failed — {e}')

n_before = adata.n_obs
adata = adata[~adata.obs['predicted_doublet']].copy()
print(f'\nAfter doublet removal:{adata.n_obs:,} cells'
      f'(removed {n_before-adata.n_obs:,})')

os.makedirs(os.path.dirname(OUT_PATH),exist_ok=True)
adata.write(OUT_PATH)
print(f'Wrote {OUT_PATH}')

#Reclaim Disk
if os.path.exists('data/interim/02_qc.h5ad'):
    os.remove('data/interim/02_qc.h5ad')
    print('Removed previous intermediate')