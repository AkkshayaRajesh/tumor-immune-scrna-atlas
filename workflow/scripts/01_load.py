"""01_load.py — Load GSE171145 into a single AnnData."""

import os
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

SAMPLE_TO_GSM = {'LJQ':   'GSM5219674', 'GBG':   'GSM5219675',
    'LYB.1': 'GSM5219676', 'LYB.2': 'GSM5219677',
    'CYD':   'GSM5219678', 'CYZ':   'GSM5219679',
    'XMS':   'GSM5219680', 'ZYQ':   'GSM5219681',
    'TGS':   'GSM5219682',
    }

RAW_DIR = 'data/raw'
META_CSV = 'data/sample_metadata.csv'

try:
    OUT_PATH = snakemake.output[0]   # noqa: F821
except NameError:
    OUT_PATH = 'data/raw/01_concatenated.h5ad'

def find_file(gsm,sample,suffix):
    expected = f'{gsm}_{sample}-T.{suffix}.gz'
    path = os.path.join(RAW_DIR, expected)
    if os.path.exists(path):
        return path
    for fn in os.listdir(RAW_DIR):
        if fn.startswith(gsm) and fn.endswith(suffix+'.gz'):
            return os.path.join(RAW_DIR,fn)
    raise FileNotFoundError(f'No file for {gsm}/{sample}/{suffix}')
    
def load_one(sample,gsm):
    counts_path = find_file(gsm,sample,'counts.tsv')
    cellname_path = find_file(gsm,sample,'cellname.list.txt')

    counts = pd.read_csv(counts_path,sep='\t',index_col=0)
    a = ad.AnnData(counts.T)
    a.X = csr_matrix(a.X.astype('float32')) 
    a.var_names_make_unique()

    cellmap = pd.read_csv(cellname_path,sep='\t')
    cellmap = cellmap.set_index('CellIndex')['CellName']
    a.obs['barcode'] = [cellmap.get(idx, idx) for idx in a.obs_names]
    a.obs['sample_id'] = sample
    a.obs_names = [f'{sample}_{bc}' for bc in a.obs['barcode']]
    return a

def main():
    adatas={}
    for sample,gsm in SAMPLE_TO_GSM.items():
        a = load_one(sample,gsm)
        adatas[sample]=a
        print(f'{sample:7s} ({gsm}): {a.n_obs:6,} cells x {a.n_vars:6,} genes')

    adatas = ad.concat(adatas, join='outer',label='sample_id_check')
    if 'sample_id_check' in adatas.obs:
        del adatas.obs['sample_id_check']
    print(f'\nConcatenated: {adatas.n_obs:,} cells × {adatas.n_vars:,} genes')

    if os.path.exists(META_CSV):
        meta = pd.read_csv(META_CSV)
        obs = adatas.obs.reset_index().merge(meta, on='sample_id', how='left')
        obs = obs.set_index('index')
        obs.index.name = None
        adatas.obs = obs
        print('\nEGFR distribution (cells):')
        print(adatas.obs['egfr_status'].value_counts().to_string())
 
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    adatas.write(OUT_PATH)
    print(f'\nWrote {OUT_PATH}')
 
 
if __name__ == '__main__':
    main()