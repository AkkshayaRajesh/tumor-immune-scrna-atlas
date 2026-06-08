"""09_de_wilcoxon.py — Per-cell-type Wilcoxon DE: EGFR mutant vs WT."""

import os
import scanpy as sc
 
adata = sc.read_h5ad('data/processed/annotated.h5ad')
os.makedirs('results/tables/de_wilcoxon', exist_ok=True)

de_results = {}
for ct in adata.obs['cell_type'].unique():
    sub = adata[adata.obs['cell_type'] == ct].copy()
    if sub.obs['egfr_status'].nunique() < 2:
        continue
    if sub.n_obs <50:
        print(f'{ct}: skipped ({sub.n_obs} cells)')
        continue
    sc.tl.rank_genes_groups(
        sub, 'egfr_status', method='wilcoxon', use_raw=True)
    df_mut = sc.get.rank_genes_groups_df(sub, group='mutant')
    df_mut.to_csv(f'results/tables/de_wilcoxon/de_{ct}.csv', index=False)
    de_results[ct] = df_mut
    n_sig = ((df_mut['pvals_adj'] < 0.05) & (df_mut['logfoldchanges'].abs() > 0.5)).sum()
    print(f'{ct}: {n_sig} genes at adj-p < 0.05 & |LFC| > 0.5')
        


#AFTER CORRECTION
'''
Biological Interpretation of EGFR Mutant vs. WT DE
Malignant Epithelial (4,893 genes) — Intrinsic oncogenic signature
This is the direct effect of constitutively active EGFR signaling: hyperactivation of PI3K/AKT/mTOR and RAS/MAPK/ERK cascades, altered cell cycle, survival, and EMT programs. Serves as your positive control — expected to be strongest.

Fibroblasts (1,809 genes, 2nd strongest) — Stromal reprogramming
EGFR-mutant tumor cells secrete TGF-β, HGF, and FGF that activate cancer-associated fibroblasts (CAFs). This large signal suggests EGFR-mutant tumors create a distinctly remodeled stroma — potentially desmoplastic or immunosuppressive — compared to WT. Clinically relevant because CAF activation can physically restrict T cell infiltration.

T_CD8 (1,615 genes) — Altered cytotoxic T cell states
EGFR-mutant tumors are known to suppress CD8 T cell recruitment and function via:

Reduced CXCL9/CXCL10 chemokine production (less homing)
Increased PD-L1 and TGF-β expression
Differential exhaustion programs
This likely explains at the transcriptional level why EGFR-mutant NSCLC patients respond poorly to checkpoint blockade.

B cells (1,544 genes, 93% robust) — Tertiary lymphoid structure differences
The high robustness here is notable — nearly all significant genes have large effect sizes. EGFR-mutant vs. WT tumors differ substantially in tertiary lymphoid structure (TLS) formation and B cell activation states. This is an emerging area worth pathway follow-up; TLS density is a positive prognostic marker and an immunotherapy response predictor.

TAM (1,225 genes) — Macrophage polarization
EGFR-mutant tumor cells drive M2-like (immunosuppressive) TAM polarization through IL-10, TGF-β, and M-CSF secretion. The drop from 2,272 → 1,225 after LFC filtering means about half the original signal was low-magnitude noise, but the remaining ~1,200 genes represent genuine macrophage reprogramming.

T_CD4 (886 genes, only 47% robust) — Interpret cautiously
More than half the T_CD4 signal was a statistical artifact of large cell numbers. The filtered 886 genes represent real helper T cell state differences — Th1/Th2/Th17 balance, CD4 support for CD8 responses — but this compartment's biology is less robustly captured in this dataset than it appeared initially.

Endothelial + Monocyte (200/287 genes, ~98% robust) — Small but genuine signals
Nearly all hits pass the LFC filter, so these are real. Endothelial differences likely reflect EGFR-driven VEGF regulation and angiogenic remodeling. Monocyte differences may reflect altered differentiation trajectories into TAMs in EGFR-mutant TME.

pDC, NK, Mast — Modest microenvironmental shaping
Smaller signals consistent with secondary effects of the EGFR-mutant TME on innate immune compartments. NK downregulation could reflect altered NKG2D ligand expression on tumor cells; pDC differences may reflect altered type I interferon tone.

MDSC (20) and Treg (2) — Almost certainly underpowered, not biologically null
Both cell types are classically enriched in immunosuppressive TMEs. Near-zero hits almost certainly reflect marginal cell counts rather than no EGFR effect. Do not interpret as "no difference" — check N before drawing conclusions.'''

# BEFORE LFC Filter/Correction
'''Results
python .\workflow\scripts\09_de_wilcoxon.py
T_CD8: 3009 genes at adj-p < 0.05
B: 1666 genes at adj-p < 0.05
Plasma: 1286 genes at adj-p < 0.05
T_CD4: 1868 genes at adj-p < 0.05
Fibroblast: 2365 genes at adj-p < 0.05
NK: 271 genes at adj-p < 0.05
Malignant_epi: 6542 genes at adj-p < 0.05
pDC: 94 genes at adj-p < 0.05
TAM: 2272 genes at adj-p < 0.05
MDSC: 21 genes at adj-p < 0.05
Mast: 79 genes at adj-p < 0.05
Endothelial: 204 genes at adj-p < 0.05
Monocyte: 292 genes at adj-p < 0.05
Treg: 2 genes at adj-p < 0.05'''

'''Takeaways
1. Malignannt Epi - Strongest
2. large signals in CD8 T cells and fibroblasts
'''