#!/usr/bin/env python
"""Compute HDAC4 repression-target activity scores two independent ways.

Method 1: decoupler ULM with a hand-curated MEF2C/RUNX2 regulon restricted to the
          repression-target set.
Method 2: scanpy.tl.score_genes scaled mean expression over the same set.

Both methods score TARGET EXPRESSION. Inferred HDAC4 repressive activity is the
negative of that (fixed sign convention, PREREGISTRATION section 4).

Outputs:
  data/processed/<tag>.scored.h5ad
  results/tables/pseudobulk_scores.tsv         (one row per sample x zone)
  results/tables/gene_set_coverage.tsv
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

DATASETS = ["human_GSE288028", "mouse_GSE288028", "mouse_GSE288529"]

SETS = {
    "primary": C.TARGETS_PRIMARY,
    "secondary": C.TARGETS_SECONDARY,
    "housekeeping": C.HOUSEKEEPING,
}


def build_regulon(genes: list[str], source: str = "MEF2C_RUNX2") -> pd.DataFrame:
    """Hand-curated regulon: HDAC4 blocks MEF2C/RUNX2, so all targets carry weight +1.

    A positive ULM score means the MEF2C/RUNX2 programme is ON, i.e. HDAC4 repressive
    activity is inferred to be LOW.
    """
    return pd.DataFrame({"source": source, "target": genes, "weight": 1.0})


def score_dataset(ad: sc.AnnData, species: str) -> tuple[sc.AnnData, list[dict]]:
    import decoupler as dc

    ad.X = ad.layers["lognorm"].copy()
    coverage = []

    for name, human_genes in SETS.items():
        mapped = C.map_genes(human_genes, species)
        found = [g for g in mapped if g in ad.var_names]
        missing = [g for g in mapped if g not in ad.var_names]
        coverage.append(
            {
                "dataset_species": species,
                "gene_set": name,
                "n_requested": len(mapped),
                "n_found": len(found),
                "found": ",".join(found),
                "missing": ",".join(missing),
            }
        )
        if not found:
            ad.obs[f"target_{name}_scanpy"] = np.nan
            ad.obs[f"target_{name}_ulm"] = np.nan
            continue

        # --- method 2: scaled mean expression
        sc.tl.score_genes(ad, found, score_name=f"target_{name}_scanpy", random_state=C.SEED)

        # --- method 1: decoupler ULM on the curated regulon
        net = build_regulon(found)
        # decoupler densifies each batch, so cap the batch so that
        # bsize x n_vars x 8 bytes stays near 1.6 GB rather than the whole matrix.
        bsize = max(1000, int(2e8 / max(ad.n_vars, 1)))
        try:
            dc.mt.ulm(ad, net, tmin=min(3, len(found)), verbose=False, bsize=bsize)
            key = next(k for k in ("score_ulm", "ulm_estimate") if k in ad.obsm)
            est = ad.obsm[key]
            col = est.columns[0] if hasattr(est, "columns") else 0
            ad.obs[f"target_{name}_ulm"] = np.asarray(
                est[col] if hasattr(est, "columns") else est[:, 0]
            )
            for k in list(ad.obsm.keys()):
                if k.startswith(("score_", "padj_", "pvals_")) and k.endswith(
                    ("ulm", "_ulm")
                ):
                    del ad.obsm[k]
        except Exception as exc:  # noqa: BLE001 - record, do not silently drop a method
            print(f"    [decoupler ULM failed for {name}: {exc}]", flush=True)
            ad.obs[f"target_{name}_ulm"] = np.nan

    # pathway components: raw expression, never folded into the score
    for g in C.map_genes(C.PATHWAY_COMPONENTS, species):
        if g in ad.var_names:
            v = ad[:, g].X
            ad.obs[f"expr_{g}"] = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.asarray(v).ravel()

    # inferred HDAC4 repressive activity = -(target expression score)
    for name in SETS:
        for meth in ("scanpy", "ulm"):
            ad.obs[f"hdac4_activity_{name}_{meth}"] = -ad.obs[f"target_{name}_{meth}"]

    return ad, coverage


def pseudobulk(ad: sc.AnnData, tag: str) -> pd.DataFrame:
    score_cols = [
        f"target_{n}_{m}" for n in SETS for m in ("scanpy", "ulm")
    ] + [f"hdac4_activity_{n}_{m}" for n in SETS for m in ("scanpy", "ulm")]
    expr_cols = [c for c in ad.obs.columns if c.startswith("expr_")]
    tech_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "S_score", "G2M_score"]

    keys = ["dataset", "sample_id", "donor", "arm", "species", "cellranger", "zone"]
    g = ad.obs.groupby(keys, observed=True)
    out = g[score_cols + expr_cols + tech_cols].mean()
    out["n_cells"] = g.size()
    out = out.reset_index()
    out.insert(0, "dataset_tag", tag)
    out["below_min_cells_flag"] = out["n_cells"] < C.MIN_CELLS_PER_ZONE
    return out


def main() -> int:
    pbs, cov = [], []
    for tag in DATASETS:
        p = C.PROC / f"{tag}.h5ad"
        print(f"\n=== {tag} ===", flush=True)
        ad = sc.read_h5ad(p)
        species = ad.obs["species"].iloc[0]
        ad, c = score_dataset(ad, species)
        for row in c:
            row["dataset_tag"] = tag
            print(f"    {row['gene_set']}: {row['n_found']}/{row['n_requested']} found"
                  + (f" | missing: {row['missing']}" if row["missing"] else ""), flush=True)
        cov.extend(c)
        pbs.append(pseudobulk(ad, tag))
        ad.write(C.PROC / f"{tag}.scored.h5ad", compression="gzip")

    pd.concat(pbs).to_csv(C.TABLES / "pseudobulk_scores.tsv", sep="\t", index=False)
    pd.DataFrame(cov).to_csv(C.TABLES / "gene_set_coverage.tsv", sep="\t", index=False)
    print("\ndone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
