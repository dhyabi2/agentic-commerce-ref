#!/usr/bin/env python
"""Phase D: genome-wide donor-level zonal expression atlas from GSE288529 (audited single library).
For ALL genes, per growth-plate STATE: mean log-norm, detection fraction, raw pseudobulk/10k, and
state-specificity. SINGLE LIBRARY / ONE ANIMAL -> descriptive only; no biological CIs (cells are not
biological replicates). Zone labels are the frozen Leiden->state map (annotated WITHOUT the query genes
driving the map beyond canonical markers)."""
import anndata, numpy as np, pandas as pd, scipy.sparse as sp, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
H5 = "/home/user/agentic-commerce/growth_plate_pth1r/data/processed.h5ad"
OUT = ROOT / "results" / "tables"; OUT.mkdir(parents=True, exist_ok=True)

# frozen Leiden -> state map (from the prior GSE288529 audit; growth-plate states only)
STATE_MAP = {"8":"Resting","4":"Proliferative","5":"Proliferative","16":"Proliferative",
             "3":"Prehypertrophic","7":"Prehypertrophic",
             "6":"Hypertrophic","12":"Hypertrophic","17":"Hypertrophic","19":"Hypertrophic"}
STATES = ["Resting","Proliferative","Prehypertrophic","Hypertrophic"]

a = anndata.read_h5ad(H5)
a.obs["state"] = a.obs["leiden"].astype(str).map(STATE_MAP)
gp = a[a.obs["state"].isin(STATES)].copy()
genes = gp.var_names.to_numpy()
ln = gp.layers["lognorm"]; cnt = gp.layers["counts"]
ln = sp.csr_matrix(ln); cnt = sp.csr_matrix(cnt)
tot = np.asarray(cnt.sum(1)).ravel()   # total counts per cell

rows_expr = []
spec = {}
per_state_mean = {}
for s in STATES:
    m = (gp.obs["state"] == s).values
    n = int(m.sum())
    ln_s = ln[m]; cnt_s = cnt[m]; tot_s = tot[m].sum()
    mean_ln = np.asarray(ln_s.mean(0)).ravel()
    det = np.asarray((ln_s > 0).sum(0)).ravel() / n
    pb = np.asarray(cnt_s.sum(0)).ravel() / max(tot_s, 1) * 1e4
    per_state_mean[s] = mean_ln
    for gi, g in enumerate(genes):
        rows_expr.append((g, s, n, float(mean_ln[gi]), float(det[gi]), float(pb[gi])))

expr = pd.DataFrame(rows_expr, columns=["gene","state","n_cells","mean_lognorm","detection_frac","pseudobulk_cp10k"])
expr["library"] = "GSE288529 GSM8769462 (single 4wk female mouse)"
expr["descriptive_only"] = True
expr.to_parquet(OUT / "GENOME_WIDE_ZONE_EXPRESSION.parquet", index=False)

# state specificity: which zone has max mean expression + a specificity score (max / sum)
M = np.vstack([per_state_mean[s] for s in STATES])          # 4 x genes
tot_mean = M.sum(0) + 1e-9
argmax = M.argmax(0)
spec_df = pd.DataFrame({
    "gene": genes,
    "max_state": [STATES[i] for i in argmax],
    "specificity_score": (M.max(0) / tot_mean),
    **{f"mean_{s}": per_state_mean[s] for s in STATES},
    "resting_vs_hypertrophic_log2fc": np.log2((per_state_mean["Resting"]+1e-3)/(per_state_mean["Hypertrophic"]+1e-3)),
})
spec_df.to_parquet(OUT / "GENOME_WIDE_STATE_SPECIFICITY.parquet", index=False)

print(f"genes: {len(genes)}; states: {STATES}")
print("wrote GENOME_WIDE_ZONE_EXPRESSION.parquet ({} rows)".format(len(expr)))
print("wrote GENOME_WIDE_STATE_SPECIFICITY.parquet ({} genes)".format(len(spec_df)))
# quick sanity: top resting-specific genes (detected)
top = spec_df[(spec_df.max_state=="Resting")].sort_values("specificity_score", ascending=False)
top = top[spec_df["mean_Resting"]>0.5].head(10)
print("\ntop resting-specific genes:", top["gene"].tolist())
