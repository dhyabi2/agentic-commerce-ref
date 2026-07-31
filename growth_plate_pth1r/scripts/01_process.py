#!/usr/bin/env python
"""
GSE288529 growth-plate scRNA-seq: QC -> normalization -> PCA -> neighbors -> UMAP -> Leiden.
Saves processed AnnData and prints per-cluster marker means for annotation.
Pth1r is deliberately NOT used to drive clustering/annotation.
"""
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 1
np.random.seed(0)

H5 = "data/GSM8769462_3402_filtered_feature_bc_matrix.h5"

adata = sc.read_10x_h5(H5)
adata.var_names_make_unique()
print("Loaded:", adata.shape)

# ---- QC metrics ----
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.match(r"^Rp[sl]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=None)
print("Pre-filter QC summary:")
print(adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].describe())

# ---- Basic filtering ----
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 15].copy()
# remove obvious high-count doublet-like outliers
upper = np.percentile(adata.obs["total_counts"], 99.5)
adata = adata[adata.obs["total_counts"] < upper].copy()
print("Post-filter:", adata.shape)

# ---- Preserve raw counts ----
adata.layers["counts"] = adata.X.copy()

# ---- Normalize + log ----
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()   # store log-normalized for quantification later

# ---- HVG on log-normalized ----
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
adata.raw = adata  # keep full log-normalized for plotting/marker lookup

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.uns["pca"] = adata_hvg.uns["pca"]

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph", n_iterations=2, directed=False)
print("Leiden clusters:", adata.obs["leiden"].value_counts().sort_index().to_dict())

# ---- Per-cluster marker means (log-normalized) for annotation ----
panels = {
    "Resting": ["Apoe", "Clu", "C1qtnf3", "Pthlh"],
    "Proliferative": ["Acan", "Col2a1", "Sox9", "Ccnd1"],
    "Prehypertrophic": ["Ihh", "Runx2", "Sp7"],
    "Hypertrophic": ["Col10a1", "Mmp13", "Ibsp"],
    "Articular": ["Prg4", "Gdf5", "Cytl1"],
    "Osteo/Perichondrial": ["Col1a1", "Col3a1", "Postn", "Bglap", "Alpl"],
    "Immune": ["Ptprc", "C1qa", "Cd68", "Lyz2"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr"],
    "Pth1r(diagnostic)": ["Pth1r"],
}
all_genes = sorted({g for v in panels.values() for g in v})
present = [g for g in all_genes if g in adata.raw.var_names]
missing = [g for g in all_genes if g not in adata.raw.var_names]
print("Marker genes missing from data:", missing)

expr = pd.DataFrame(
    adata.raw[:, present].X.toarray(),
    columns=present,
    index=adata.obs_names,
)
expr["leiden"] = adata.obs["leiden"].values
cluster_means = expr.groupby("leiden", observed=True).mean()
pd.set_option("display.width", 220, "display.max_columns", 60)
print("\n=== Per-cluster mean log-normalized expression ===")
print(cluster_means.round(2).to_string())

# percent detected per cluster for key markers
det = (expr.drop(columns="leiden") > 0).astype(int)
det["leiden"] = adata.obs["leiden"].values
cluster_pct = det.groupby("leiden", observed=True).mean() * 100
print("\n=== Per-cluster percent detected (%) ===")
print(cluster_pct.round(1).to_string())

cluster_means.to_csv("tables/cluster_marker_means.csv")
cluster_pct.to_csv("tables/cluster_marker_pctdetected.csv")
adata.write("data/processed.h5ad")
print("\nSaved data/processed.h5ad")
