import scanpy as sc
samples = {}
for sample_id in ['LJQ', 'GBG', 'LYB.1', 'LYB.2', 'CYD', 'CYZ', 'XMS', 'ZYQ', 'TGS']:
    samples[sample_id] = sc.read_10x_mtx(f'data/raw/{sample_id}/')
adata = sc.concat(samples, label='patient_id')