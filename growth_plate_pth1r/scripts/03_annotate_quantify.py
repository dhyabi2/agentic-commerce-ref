#!/usr/bin/env python
"""
Annotate growth-plate states (Pth1r-independent), quantify Pth1r per state,
compare Resting vs Prehypertrophic, test for a Pth1r+ resting subpopulation,
and generate figures + CSV summary.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import ranksums, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
sc.settings.figdir = "figures"
sc.settings.verbosity = 0

adata = sc.read_h5ad("data/processed.h5ad")
adata.X = adata.layers["lognorm"].copy()          # ensure X = log-normalized
adata.raw = adata

# ---------------------------------------------------------------------------
# Final cluster -> label map. Assignment uses the marker panels ONLY
# (Pth1r never used). Rationale recorded in REPORT.md.
# ---------------------------------------------------------------------------
STATE_MAP = {
    # growth-plate states
    "8": "Resting",
    "4": "Proliferative", "5": "Proliferative", "16": "Proliferative",
    "3": "Prehypertrophic", "7": "Prehypertrophic",
    "6": "Hypertrophic", "12": "Hypertrophic", "17": "Hypertrophic", "19": "Hypertrophic",
    # excluded contaminants
    "10": "Excl:Articular", "9": "Excl:Articular", "15": "Excl:Articular",
    "0": "Excl:Osteo/Perichondrial", "14": "Excl:Osteo/Perichondrial",
    "1": "Excl:Endothelial",
    "13": "Excl:Immune",
    # low-quality / mixed / border -> unassigned, excluded from state stats
    "2": "Excl:LowQual/Mixed", "11": "Excl:LowQual/Mixed", "18": "Excl:LowQual/Mixed",
    "20": "Excl:LowQual/Mixed", "21": "Excl:LowQual/Mixed", "22": "Excl:LowQual/Mixed",
    "23": "Excl:LowQual/Mixed",
}
adata.obs["state"] = adata.obs["leiden"].astype(str).map(STATE_MAP).astype("category")
print("State composition (all cells):")
print(adata.obs["state"].value_counts().to_string())

STATES = ["Resting", "Proliferative", "Prehypertrophic", "Hypertrophic"]

def gene_vec(ad, gene, layer=None):
    X = ad.layers[layer] if layer else ad.X
    j = ad.var_names.get_loc(gene)
    v = X[:, j]
    return np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()

# ---------------------------------------------------------------------------
# 6. Pth1r quantification per annotated state
# ---------------------------------------------------------------------------
ln = gene_vec(adata, "Pth1r", "lognorm")
raw = gene_vec(adata, "Pth1r", "counts")
total = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()

rows = []
for s in STATES:
    m = (adata.obs["state"] == s).values
    ln_s, raw_s, tot_s = ln[m], raw[m], total[m]
    pseudobulk_cp10k = raw_s.sum() / tot_s.sum() * 1e4
    rows.append({
        "state": s,
        "n_cells": int(m.sum()),
        "mean_lognorm": float(ln_s.mean()),
        "median_lognorm": float(np.median(ln_s)),
        "pct_detected": float((ln_s > 0).mean() * 100),
        "pseudobulk_raw_per_10k": float(pseudobulk_cp10k),
    })
summary = pd.DataFrame(rows).set_index("state")
pd.set_option("display.width", 200)
print("\n=== Pth1r per state ===")
print(summary.round(4).to_string())
summary.to_csv("tables/pth1r_by_state.csv")

# ---------------------------------------------------------------------------
# 7. Resting vs Prehypertrophic
# ---------------------------------------------------------------------------
rest_m = (adata.obs["state"] == "Resting").values
pre_m  = (adata.obs["state"] == "Prehypertrophic").values
ln_rest, ln_pre = ln[rest_m], ln[pre_m]
raw_rest, raw_pre = raw[rest_m], raw[pre_m]
tot_rest, tot_pre = total[rest_m], total[pre_m]

expr_ratio = ln_pre.mean() / ln_rest.mean()
det_rest = (ln_rest > 0).mean()
det_pre  = (ln_pre > 0).mean()
det_ratio = det_pre / det_rest
pb_rest = raw_rest.sum() / tot_rest.sum() * 1e4
pb_pre  = raw_pre.sum() / tot_pre.sum() * 1e4

w_stat, w_p = ranksums(ln_pre, ln_rest)          # Wilcoxon rank-sum
u_stat, u_p = mannwhitneyu(ln_pre, ln_rest, alternative="two-sided")

def boot_ci(a, b, fn, n=5000):
    vals = np.empty(n)
    na, nb = len(a), len(b)
    for i in range(n):
        vals[i] = fn(a[rng.integers(0, na, na)], b[rng.integers(0, nb, nb)])
    return np.percentile(vals, [2.5, 97.5])

ci_expr = boot_ci(ln_pre, ln_rest, lambda p, r: p.mean() / max(r.mean(), 1e-9))
ci_det  = boot_ci(ln_pre, ln_rest, lambda p, r: (p > 0).mean() / max((r > 0).mean(), 1e-9))
ci_diff = boot_ci(ln_pre, ln_rest, lambda p, r: p.mean() - r.mean())

comp = {
    "resting_n": int(rest_m.sum()), "prehyp_n": int(pre_m.sum()),
    "resting_mean_lognorm": float(ln_rest.mean()), "prehyp_mean_lognorm": float(ln_pre.mean()),
    "expr_ratio_pre_over_rest": float(expr_ratio),
    "expr_ratio_95CI_low": float(ci_expr[0]), "expr_ratio_95CI_high": float(ci_expr[1]),
    "mean_diff_pre_minus_rest": float(ln_pre.mean() - ln_rest.mean()),
    "mean_diff_95CI_low": float(ci_diff[0]), "mean_diff_95CI_high": float(ci_diff[1]),
    "resting_pct_detected": float(det_rest * 100), "prehyp_pct_detected": float(det_pre * 100),
    "detection_ratio_pre_over_rest": float(det_ratio),
    "detection_ratio_95CI_low": float(ci_det[0]), "detection_ratio_95CI_high": float(ci_det[1]),
    "resting_pseudobulk_cp10k": float(pb_rest), "prehyp_pseudobulk_cp10k": float(pb_pre),
    "pseudobulk_ratio_pre_over_rest": float(pb_pre / pb_rest),
    "wilcoxon_ranksum_stat": float(w_stat), "wilcoxon_ranksum_p": float(w_p),
    "mannwhitneyu_stat": float(u_stat), "mannwhitneyu_p": float(u_p),
}
print("\n=== Resting vs Prehypertrophic ===")
for k, v in comp.items():
    print(f"{k:38s}: {v:.4g}" if isinstance(v, float) else f"{k:38s}: {v}")
pd.Series(comp).to_csv("tables/resting_vs_prehyp.csv", header=["value"])

# ---------------------------------------------------------------------------
# 8. Is Pth1r+ within Resting a distinct subpopulation, or randomly scattered?
#    (a) KNN homophily: fraction of a Pth1r+ resting cell's neighbours (within
#        resting) that are also Pth1r+, vs a label-permutation null.
#    (b) Moran's I of Pth1r on the resting-only kNN graph.
# ---------------------------------------------------------------------------
sub = adata[rest_m].copy()
sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca", random_state=0)
pos = (gene_vec(sub, "Pth1r", "lognorm") > 0).astype(float)
A = sub.obsp["connectivities"]
A = (A > 0).astype(float)                    # binarize adjacency
deg = np.asarray(A.sum(1)).ravel()
deg[deg == 0] = 1
neigh_posfrac = np.asarray(A.dot(pos)).ravel() / deg
obs_homophily = neigh_posfrac[pos == 1].mean()   # among Pth1r+ cells

null = np.empty(2000)
idx = np.arange(len(pos))
for i in range(2000):
    perm = pos[rng.permutation(idx)]
    nf = np.asarray(A.dot(perm)).ravel() / deg
    null[i] = nf[perm == 1].mean()
null_mean, null_sd = null.mean(), null.std()
z_homo = (obs_homophily - null_mean) / (null_sd + 1e-12)
p_homo = (np.sum(null >= obs_homophily) + 1) / (len(null) + 1)

# Moran's I of continuous Pth1r within resting
x = gene_vec(sub, "Pth1r", "lognorm")
xc = x - x.mean()
W = A
Wsum = W.sum()
num = xc.dot(W.dot(xc))
den = (xc ** 2).sum()
morans_I = (len(x) / Wsum) * (num / den)
null_I = np.empty(2000)
for i in range(2000):
    xp = xc[rng.permutation(len(xc))]
    null_I[i] = (len(x) / Wsum) * (xp.dot(W.dot(xp)) / den)
p_I = (np.sum(null_I >= morans_I) + 1) / (len(null_I) + 1)

subpop = {
    "resting_pth1r_pos_frac": float(pos.mean()),
    "observed_knn_homophily": float(obs_homophily),
    "null_homophily_mean": float(null_mean), "null_homophily_sd": float(null_sd),
    "homophily_z": float(z_homo), "homophily_perm_p": float(p_homo),
    "morans_I_pth1r_resting": float(morans_I),
    "morans_I_null_mean": float(null_I.mean()), "morans_I_perm_p": float(p_I),
}
print("\n=== Pth1r+ resting subpopulation test ===")
for k, v in subpop.items():
    print(f"{k:32s}: {v:.4g}")
pd.Series(subpop).to_csv("tables/resting_subpopulation.csv", header=["value"])

# ---------------------------------------------------------------------------
# 9. Figures
# ---------------------------------------------------------------------------
adata.obs["state_plot"] = adata.obs["state"].astype(str)
sc.pl.umap(adata, color="state_plot", title="Annotated states", show=False,
           save="_states.png", size=12, legend_fontsize=7)
sc.pl.umap(adata, color="leiden", title="Leiden clusters", show=False,
           save="_leiden.png", size=12, legend_loc="on data", legend_fontsize=6)
for g in ["Apoe", "Pthlh", "Ihh", "Col10a1", "Pth1r"]:
    sc.pl.umap(adata, color=g, title=g, show=False, save=f"_{g}.png", size=12, cmap="viridis")

# dot plot markers by cluster
dot_markers = ["Apoe", "Clu", "C1qtnf3", "Pthlh", "Acan", "Col2a1", "Sox9", "Ccnd1",
               "Ihh", "Runx2", "Sp7", "Col10a1", "Mmp13", "Ibsp", "Pth1r"]
sc.pl.dotplot(adata, dot_markers, groupby="leiden", show=False, save="markers_by_cluster.png",
              standard_scale="var")

# violin Pth1r by state (growth-plate states only)
gp = adata[adata.obs["state"].isin(STATES)].copy()
gp.obs["state"] = pd.Categorical(gp.obs["state"].astype(str), categories=STATES, ordered=True)
sc.pl.violin(gp, "Pth1r", groupby="state", rotation=30, show=False, save="_Pth1r_by_state.png")

# dotplot Pth1r by state
sc.pl.dotplot(gp, ["Pth1r"], groupby="state", show=False, save="Pth1r_by_state.png")

print("\nFigures written to figures/. Done.")
