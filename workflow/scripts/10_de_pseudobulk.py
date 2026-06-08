"""10_de_pseudobulk.py — Patient-level pseudobulk DE with DESeq2."""
import os
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

adata = sc.read_h5ad('data/processed/annotated.h5ad')
os.makedirs('results/tables/de_pseudobulk', exist_ok=True)

# Use full gene set from .raw for DE
raw = adata.raw.to_adata()
raw.obs = adata.obs.copy()

for ct in adata.obs['cell_type'].unique():
    sub = raw[raw.obs['cell_type']==ct].copy()
    if sub.obs['egfr_status'].nunique() < 2 or sub.n_obs <50:
        continue
    
    # Pseudobulk: sum counts per (sample_id, cell_type)
    # .raw.to_adata() never carries layers, so .X is always the raw count matrix
    mat = sub.X
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    pb = (
        pd.DataFrame(mat, index=sub.obs_names, columns=sub.var_names)
        .groupby(sub.obs['sample_id'].values).sum()
    )
    metadata = sub.obs.drop_duplicates('sample_id').set_index(
        'sample_id')[['egfr_status']].loc[pb.index]
    
    if metadata['egfr_status'].nunique() < 2: continue
    # DESeq2 needs ≥3 replicates per group to estimate within-group dispersion;
    # 2-per-group leaves 0 degrees of freedom and causes convergence failures.
    if (metadata['egfr_status'].value_counts() < 3).any(): continue

    # Drop genes with < 10 total pseudobulk counts — too sparse for reliable DE
    pb = pb.loc[:, pb.sum(axis=0) >= 10]
    if pb.shape[1] == 0:
        print(f'{ct}: skipped — no genes passed count filter')
        continue

    try:
        dds = DeseqDataSet(
            counts=pb.astype(int),
            metadata=metadata,
            design='~ egfr_status',
            refit_cooks=True,
            quiet=True,
        )
        dds.deseq2()
        ds = DeseqStats(dds, contrast=['egfr_status', 'mutant', 'wildtype'],
                        quiet=True)
        ds.summary()
        res = ds.results_df.sort_values('padj').reset_index()
        res.to_csv(f'results/tables/de_pseudobulk/de_{ct}.csv', index=False)
        n_sig = (res['padj'] < 0.05).sum()
        counts = metadata['egfr_status'].value_counts().to_dict()
        print(f'{ct}: {n_sig} genes at padj < 0.05 (pseudobulk) | n_samples={counts}')
    except Exception as e:
        print(f'{ct}: pseudobulk failed — {e}')

        