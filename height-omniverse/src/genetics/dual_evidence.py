#!/usr/bin/env python
"""Phase G: dual-/multi-evidence discovery table. For each gene, combine INDEPENDENT evidence modes:
 (1) human genetics direction, (2) growth-plate expression, (3) zone specificity, (4) compound exists.
Component scores stay VISIBLE. 'high_priority_research' requires >=3 modes. NO recommendation."""
import pandas as pd, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tables"
spec = pd.read_parquet(OUT/"GENOME_WIDE_STATE_SPECIFICITY.parquet")
gd = pd.read_csv(ROOT/"data/curated/genetic_direction.csv")
STATES=["Resting","Proliferative","Prehypertrophic","Hypertrophic"]
spec["gp_expr"]=spec[[f"mean_{s}" for s in STATES]].sum(1)
def mouse(s): return s[:1].upper()+s[1:].lower()

# compound availability (Phase F) if present
comp_path = OUT/"COMPOUND_TARGET.parquet"
comp_targets = set()
if comp_path.exists():
    comp_targets = set(pd.read_parquet(comp_path)["target_gene"].str.upper())

rows=[]
for g in sorted(set(gd["gene"])):
    grow = gd[gd["gene"]==g].iloc[0]
    m = mouse(g)
    srow = spec[spec["gene"]==m]
    gp_expr = float(srow["gp_expr"].iloc[0]) if len(srow) else 0.0
    zone = srow["max_state"].iloc[0] if len(srow) else "NA(absent in mouse)"
    zspec = float(srow["specificity_score"].iloc[0]) if len(srow) else 0.0
    m_genetics = 1
    m_expr = int(gp_expr > 0.5)
    m_zone = int(zspec > 0.4 and gp_expr > 0.5)
    m_comp = int(g.upper() in comp_targets)
    modes = m_genetics + m_expr + m_zone + m_comp
    rows.append({"gene":g, "height_direction":grow["height_direction"], "pathway":grow["pathway"],
                 "gp_expression":round(gp_expr,2), "top_zone":zone, "zone_specificity":round(zspec,2),
                 "mode_genetics":m_genetics, "mode_gp_expressed":m_expr, "mode_zone_specific":m_zone,
                 "mode_has_compound":m_comp, "n_evidence_modes":modes,
                 "high_priority_research": modes>=3, "source":grow["source_id"]})
df = pd.DataFrame(rows).sort_values(["n_evidence_modes","gp_expression"], ascending=False)
df.to_csv(OUT/"DUAL_EVIDENCE_CANDIDATES.csv", index=False)
print(f"genes scored: {len(df)}; compound layer loaded: {bool(comp_targets)}")
print(df[["gene","height_direction","top_zone","gp_expression","n_evidence_modes","high_priority_research"]].to_string(index=False))

# negative-space tables
neg = {
 "genetics_without_gp_expression": df[df.mode_gp_expressed==0]["gene"].tolist(),
 "high_gp_expr_no_compound": df[(df.mode_gp_expressed==1)&(df.mode_has_compound==0)]["gene"].tolist(),
}
pd.Series({k:";".join(v) for k,v in neg.items()}).to_csv(OUT/"NEGATIVE_SPACE.csv", header=["genes"])
print("\nnegative space written (genetics-without-GP-expression; GP-expressed-without-compound)")
