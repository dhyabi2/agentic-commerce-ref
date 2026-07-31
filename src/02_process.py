#!/usr/bin/env python
"""Load, QC, cluster and zone-annotate the growth plate scRNA-seq datasets.

Outputs:
  data/processed/human_GSE288028.h5ad
  data/processed/mouse_GSE288028.h5ad
  data/processed/mouse_GSE288529.h5ad
  results/tables/qc_per_sample.tsv
  results/tables/cluster_zone_assignment_<dataset>.tsv
  results/tables/cells_per_zone_per_sample.tsv
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hdac4_config as C  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.n_jobs = 4
sc.settings.verbosity = 1

CHONDRO = ["COL2A1", "ACAN", "SOX9", "COL9A1", "COL11A1", "COL11A2"]
IMMUNE = ["PTPRC", "CD68", "LYZ", "CD3E"]
ENDOTHELIAL = ["PECAM1", "CDH5", "KDR"]
ERYTHROID = ["HBB", "HBA1", "ALAS2"]


def present(genes: list[str], adata) -> list[str]:
    return [g for g in genes if g in adata.var_names]


def regress_out_fast(X: np.ndarray, covars: np.ndarray) -> np.ndarray:
    """OLS-residualise every column of X against covars (plus intercept).

    Mathematically identical to scanpy.pp.regress_out but solved as one least-squares
    problem instead of a per-gene loop.
    """
    D = np.column_stack([np.ones(covars.shape[0], dtype=np.float64), covars.astype(np.float64)])
    beta, *_ = np.linalg.lstsq(D, X.astype(np.float64), rcond=None)
    return (X.astype(np.float64) - D @ beta).astype(np.float32)


def load_gse288028(species: str) -> sc.AnnData:
    rows = [r for r in C.GSE288028_SAMPLES if r[4] == species]
    parts = []
    for gsm, fname, donor, arm, sp, crv in rows:
        a = sc.read_10x_h5(C.RAW / "GSE288028" / fname)
        a.var_names_make_unique()
        a.obs["gsm"] = gsm
        a.obs["sample_id"] = f"{gsm}_{donor}_{arm}"
        a.obs["donor"] = donor
        a.obs["arm"] = arm
        a.obs["species"] = sp
        a.obs["cellranger"] = crv
        a.obs["dataset"] = "GSE288028"
        parts.append(a)
    ad = sc.concat(parts, index_unique="-", label=None)
    ad.obs_names_make_unique()
    return ad


def load_gse288529() -> sc.AnnData:
    a = sc.read_10x_h5(C.RAW / "GSE288529" / "GSM8769462_3402_filtered_feature_bc_matrix.h5")
    a.var_names_make_unique()
    a.obs["gsm"] = "GSM8769462"
    a.obs["sample_id"] = "GSM8769462_M3_4wk"
    a.obs["donor"] = "M3_4wk"
    a.obs["arm"] = "uncultured"
    a.obs["species"] = "mouse"
    a.obs["cellranger"] = "unknown"
    a.obs["dataset"] = "GSE288529"
    return a


def qc(ad: sc.AnnData, species: str, tag: str) -> tuple[sc.AnnData, pd.DataFrame]:
    mt_prefix = "MT-" if species == "human" else "mt-"
    ad.var["mt"] = ad.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

    pre = (
        ad.obs.groupby("sample_id", observed=True)
        .agg(
            n_cells_pre=("total_counts", "size"),
            median_counts_pre=("total_counts", "median"),
            median_genes_pre=("n_genes_by_counts", "median"),
            median_pct_mt_pre=("pct_counts_mt", "median"),
        )
        .reset_index()
    )

    keep = (
        (ad.obs["n_genes_by_counts"] >= C.MIN_GENES_PER_CELL)
        & (ad.obs["total_counts"] >= C.MIN_COUNTS_PER_CELL)
        & (ad.obs["pct_counts_mt"] < C.MAX_PCT_MT)
    )
    ad = ad[keep].copy()
    sc.pp.filter_genes(ad, min_cells=C.MIN_CELLS_PER_GENE)

    post = (
        ad.obs.groupby("sample_id", observed=True)
        .agg(
            n_cells_post=("total_counts", "size"),
            median_counts_post=("total_counts", "median"),
            median_genes_post=("n_genes_by_counts", "median"),
            median_pct_mt_post=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    stats = pre.merge(post, on="sample_id", how="outer")
    stats.insert(0, "dataset_tag", tag)
    return ad, stats


def cluster(ad: sc.AnnData, species: str, resolution: float = 1.0) -> sc.AnnData:
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers["lognorm"] = ad.X.copy()  # zone + activity scoring use THIS layer

    # cell cycle and stress scores (computed on unregressed lognorm)
    s = present(C.map_genes(C.S_GENES, species), ad)
    g2m = present(C.map_genes(C.G2M_GENES, species), ad)
    if s and g2m:
        sc.tl.score_genes_cell_cycle(ad, s_genes=s, g2m_genes=g2m)
    else:
        ad.obs["S_score"] = 0.0
        ad.obs["G2M_score"] = 0.0
    stress = present(C.map_genes(C.STRESS_GENES, species), ad)
    sc.tl.score_genes(ad, stress, score_name="stress_score", random_state=C.SEED)

    n_batches = ad.obs["sample_id"].nunique()
    sc.pp.highly_variable_genes(
        ad, n_top_genes=2000, batch_key="sample_id" if n_batches > 1 else None
    )
    hv = ad[:, ad.var.highly_variable].copy()

    X = hv.X.toarray() if hasattr(hv.X, "toarray") else np.asarray(hv.X)
    covars = np.column_stack(
        [
            ad.obs["S_score"].to_numpy(),
            ad.obs["G2M_score"].to_numpy(),
            ad.obs["stress_score"].to_numpy(),
            ad.obs["pct_counts_mt"].to_numpy(),
        ]
    )
    X = regress_out_fast(X, covars)
    X -= X.mean(0)
    sd = X.std(0)
    sd[sd == 0] = 1.0
    X /= sd
    np.clip(X, -10, 10, out=X)
    hv.X = X

    sc.tl.pca(hv, n_comps=50, svd_solver="arpack", random_state=C.SEED)
    sc.pp.neighbors(hv, n_neighbors=15, n_pcs=50, random_state=C.SEED)
    sc.tl.leiden(hv, resolution=resolution, key_added="leiden", random_state=C.SEED, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(hv, random_state=C.SEED)

    ad.obs["leiden"] = hv.obs["leiden"].values
    ad.obsm["X_pca"] = hv.obsm["X_pca"]
    ad.obsm["X_umap"] = hv.obsm["X_umap"]
    ad.uns["regressed_covariates"] = ["S_score", "G2M_score", "stress_score", "pct_counts_mt"]
    return ad


def assign_zones(ad: sc.AnnData, species: str, tag: str) -> tuple[sc.AnnData, pd.DataFrame]:
    """Assign each leiden cluster to a zone using ONLY the zone markers.

    Scoring is done on the unregressed lognorm layer so that the proliferative zone is
    still identifiable after cell-cycle regression of the embedding.
    """
    ad.X = ad.layers["lognorm"].copy()

    used = {}
    for zone, genes in C.ZONE_MARKERS.items():
        g = present(C.map_genes(genes, species), ad)
        used[zone] = g
        sc.tl.score_genes(ad, g, score_name=f"zscore_{zone}", random_state=C.SEED)
    for name, genes in [
        ("chondro", CHONDRO), ("immune", IMMUNE), ("endo", ENDOTHELIAL), ("ery", ERYTHROID)
    ]:
        g = present(C.map_genes(genes, species), ad)
        used[name] = g
        sc.tl.score_genes(ad, g, score_name=f"lineage_{name}", random_state=C.SEED)

    cols = [f"zscore_{z}" for z in C.ZONE_ORDER]
    per_cluster = ad.obs.groupby("leiden", observed=True)[
        cols + ["lineage_chondro", "lineage_immune", "lineage_endo", "lineage_ery"]
    ].mean()

    # z-score each zone column across clusters so zones with different dynamic range
    # compete fairly
    zc = (per_cluster[cols] - per_cluster[cols].mean()) / per_cluster[cols].std(ddof=0).replace(0, 1)
    per_cluster["zone"] = [C.ZONE_ORDER[i] for i in zc.to_numpy().argmax(1)]

    non_chondro = (per_cluster["lineage_chondro"] <= 0) | (
        per_cluster[["lineage_immune", "lineage_endo", "lineage_ery"]].max(1)
        > per_cluster["lineage_chondro"]
    )
    per_cluster.loc[non_chondro, "zone"] = "non_chondrocyte"
    per_cluster["n_cells"] = ad.obs.groupby("leiden", observed=True).size()

    ad.obs["zone"] = ad.obs["leiden"].map(per_cluster["zone"]).astype(str)
    ad.uns["zone_marker_genes_used"] = {k: list(v) for k, v in used.items()}

    out = per_cluster.reset_index()
    out.insert(0, "dataset_tag", tag)
    return ad, out


def main() -> int:
    qc_frames, zone_frames, counts = [], [], []

    jobs = [
        ("human_GSE288028", lambda: load_gse288028("human"), "human"),
        ("mouse_GSE288028", lambda: load_gse288028("mouse"), "mouse"),
        ("mouse_GSE288529", load_gse288529, "mouse"),
    ]

    for tag, loader, species in jobs:
        print(f"\n{'='*70}\n{tag}\n{'='*70}", flush=True)
        ad = loader()
        print(f"  loaded {ad.n_obs} cells x {ad.n_vars} genes", flush=True)
        ad, stats = qc(ad, species, tag)
        qc_frames.append(stats)
        print(f"  post-QC {ad.n_obs} cells x {ad.n_vars} genes", flush=True)
        ad = cluster(ad, species)
        print(f"  {ad.obs['leiden'].nunique()} clusters", flush=True)
        ad, zdf = assign_zones(ad, species, tag)
        zone_frames.append(zdf)
        print(ad.obs["zone"].value_counts().to_string(), flush=True)

        cnt = (
            ad.obs.groupby(["sample_id", "donor", "arm", "zone"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        cnt.insert(0, "dataset_tag", tag)
        counts.append(cnt)

        ad.write(C.PROC / f"{tag}.h5ad", compression="gzip")
        print(f"  wrote {C.PROC / (tag + '.h5ad')}", flush=True)

    pd.concat(qc_frames).to_csv(C.TABLES / "qc_per_sample.tsv", sep="\t", index=False)
    pd.concat(zone_frames).to_csv(
        C.TABLES / "cluster_zone_assignment.tsv", sep="\t", index=False
    )
    allc = pd.concat(counts)
    allc["below_min_cells_flag"] = allc["n_cells"] < C.MIN_CELLS_PER_ZONE
    allc.to_csv(C.TABLES / "cells_per_zone_per_sample.tsv", sep="\t", index=False)
    print("\ndone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
