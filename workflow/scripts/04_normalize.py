"""04_normalize.py — Normalize, log-transform, select HVGs."""
import os
import scanpy as sc
 
try:
    IN_PATH  = snakemake.input[0]    # noqa
    OUT_PATH = snakemake.output[0]   # noqa
except NameError:
    IN_PATH  = 'data/interim/03_doublet.h5ad'
    OUT_PATH = 'data/interim/04_normalized.h5ad'

adata = sc.read_h5ad(IN_PATH)
print(f'Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes')

#1. preserve the raw counts (REQ for scVI later)
adata.layers['counts'] = adata.X.copy()

#2. Normalize to 10k counts per cell
sc.pp.normalize_total(adata,target_sum=1e4)

#3. log1p Transform
sc.pp.log1p(adata)

# 4. Snapshot the full log-normalized object as .raw
#    This preserves ALL genes (not just HVGs) for DE and marker plots later
adata.raw = adata

# 5. Select HVGs — note batch_key uses sample_id, not egfr_status
#    We want genes variable within samples (technical+biological), not
#    genes that differ between EGFR groups (we want to KEEP that signal)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor='seurat_v3',
    layer='counts',
    batch_key='sample_id',)
print(f'HVGs selected: {adata.var["highly_variable"].sum():,}')

# 6. Subset to HVGs (full set still in adata.raw)
adata = adata[:, adata.var['highly_variable']].copy()

os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
adata.write(OUT_PATH)
print(f'Wrote {OUT_PATH}: {adata.n_obs:,} × {adata.n_vars:,}')
 
if os.path.exists('data/interim/03_doublet.h5ad'):
    os.remove('data/interim/03_doublet.h5ad')