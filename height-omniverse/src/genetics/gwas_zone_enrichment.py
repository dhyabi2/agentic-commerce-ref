#!/usr/bin/env python
"""Phase E: is the verified human height-gene set enriched in specific growth-plate zones?
Tests the growth-plate-enrichment claim directly, on the GSE288529 genome-wide atlas (Phase D),
via a label-permutation null. Uses the curated verified monogenic/allelic-series height genes
(mouse orthologs); full GWAS-summary-statistic ingestion is logged as partially_complete (E1).
SINGLE library -> descriptive enrichment, not a biological-replicate test."""
import pandas as pd, numpy as np, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
rng = np.random.default_rng(20260731)
OUT = ROOT / "results" / "tables"

spec = pd.read_parquet(OUT / "GENOME_WIDE_STATE_SPECIFICITY.parquet")
atlas_genes = set(spec["gene"])
gd = pd.read_csv(ROOT / "data" / "curated" / "genetic_direction.csv")

# human -> mouse ortholog (title-case); keep only those present in the mouse atlas
def mouse(sym): return sym[:1].upper() + sym[1:].lower()
height_genes_h = sorted(set(gd["gene"]))
mapped = {g: mouse(g) for g in height_genes_h}
present = {h: m for h, m in mapped.items() if m in atlas_genes}
absent = [h for h in height_genes_h if h not in present]
hg = sorted(set(present.values()))
print(f"height genes: {len(height_genes_h)} human; {len(hg)} with mouse ortholog in atlas; absent: {absent}")

STATES = ["Resting","Proliferative","Prehypertrophic","Hypertrophic"]
sub = spec[spec["gene"].isin(hg)]
# also require the gene to be actually expressed (detected) to be zone-assignable
expr_mask = sub[[f"mean_{s}" for s in STATES]].sum(1) > 0.1
sub = sub[expr_mask]
obs_counts = sub["max_state"].value_counts().reindex(STATES).fillna(0).astype(int)

# permutation null: random gene sets of equal size drawn from EXPRESSED atlas genes
expressed = spec[spec[[f"mean_{s}" for s in STATES]].sum(1) > 0.1]
pool = expressed["gene"].to_numpy()
k = len(sub)
NPERM = 5000
null = {s: np.empty(NPERM) for s in STATES}
state_of = dict(zip(expressed["gene"], expressed["max_state"]))
pool_states = np.array([state_of[g] for g in pool])
for i in range(NPERM):
    samp = pool_states[rng.integers(0, len(pool), k)]
    for s in STATES:
        null[s][i] = (samp == s).sum()

rows = []
for s in STATES:
    obs = int(obs_counts[s]); nm = null[s].mean(); nsd = null[s].std()+1e-9
    p_hi = (np.sum(null[s] >= obs)+1)/(NPERM+1)
    rows.append({"zone": s, "height_genes_assigned": obs, "expected_null_mean": round(nm,2),
                 "enrichment_ratio": round(obs/max(nm,1e-9),2), "z": round((obs-nm)/nsd,2),
                 "perm_p_enriched": round(p_hi,4)})
res = pd.DataFrame(rows)
res.to_csv(OUT / "GWAS_ZONE_ENRICHMENT.csv", index=False)
# also: are height genes more growth-plate-EXPRESSED than random? (mean total expression)
hg_expr = expressed[expressed["gene"].isin(hg)][[f"mean_{s}" for s in STATES]].sum(1).mean()
null_expr = np.array([expressed[[f"mean_{s}" for s in STATES]].sum(1).iloc[rng.integers(0,len(expressed),k)].mean() for _ in range(2000)])
p_expr = (np.sum(null_expr >= hg_expr)+1)/(2001)
print(res.to_string(index=False))
print(f"\nheight-gene mean GP expression {hg_expr:.2f} vs null {null_expr.mean():.2f}; perm_p={p_expr:.4f}")
print("assigned height genes by zone:")
for s in STATES:
    print(f"  {s}: {sorted(sub[sub.max_state==s]['gene'].tolist())}")
pd.Series({"height_gene_mean_GP_expression":hg_expr,"null_mean":float(null_expr.mean()),
           "perm_p_more_expressed":p_expr,"n_height_genes_tested":k,"absent_in_mouse":";".join(absent)}
          ).to_csv(OUT/"GWAS_ZONE_ENRICHMENT_summary.csv", header=["value"])
