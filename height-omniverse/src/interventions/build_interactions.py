#!/usr/bin/env python
"""Phase K: pairwise intervention interaction matrix by pathway convergence.
Classifies each pair as: same-target-redundant / shared-downstream / distinct / opposing.
A BLANK/'untested' cell means no direct combination evidence exists - NOT 'safe' or 'additive'.
NO stack recommendation."""
import duckdb, pandas as pd, itertools, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
con = duckdb.connect(str(ROOT/"database/height_omniverse.duckdb"), read_only=True)
ivs = con.execute("SELECT name,target,modality FROM interventions").fetchdf()
con.close()

# pathway membership (curated from targets); an intervention can touch several nodes
PATHWAY = {
 "Recombinant GH":{"GH_IGF"}, "Long-acting GH":{"GH_IGF"}, "Mecasermin (rhIGF-1)":{"GH_IGF"},
 "Aromatase inhibitor":{"estrogen"}, "Estrogen (low-dose)":{"estrogen"},
 "GnRH agonist":{"puberty_axis"}, "Oxandrolone":{"androgen"},
 "Vosoritide":{"CNP_NPR2"}, "Navepegritide (TransCon CNP)":{"CNP_NPR2"},
 "Infigratinib":{"FGFR3"}, "TYRA-300":{"FGFR3"},
 "Levothyroxine":{"thyroid"}, "Nutrition/energy supplementation":{"substrate_IGF"},
 "Abaloparatide":{"PTH1R"},
}
# downstream convergence: pathways sharing a node
CONVERGE = {("CNP_NPR2","FGFR3"):"shared-downstream (MAPK braking; CNP inhibits FGFR3/MAPK)",
            ("estrogen","GH_IGF"):"shared-downstream (estrogen modulates GH/IGF + fusion)",
            ("puberty_axis","estrogen"):"shared-downstream (GnRHa lowers sex steroids -> delays fusion)",
            ("GH_IGF","substrate_IGF"):"same-axis (IGF availability)",
            ("PTH1R","CNP_NPR2"):"distinct (both delay hypertrophy, different receptors)"}
def rel(a,b):
    pa,pb=PATHWAY.get(a,set()),PATHWAY.get(b,set())
    if pa & pb: return "same-target/pathway-redundant"
    for (x,y),lab in CONVERGE.items():
        if (x in pa and y in pb) or (y in pa and x in pb): return lab
    return "untested/unknown"

rows=[]
for a,b in itertools.combinations(sorted(ivs["name"]),2):
    r=rel(a,b)
    rows.append({"intervention_a":a,"intervention_b":b,"pathway_a":";".join(PATHWAY.get(a,{"?"})),
                 "pathway_b":";".join(PATHWAY.get(b,{"?"})),"relation":r,
                 "combination_evidence":"none_in_corpus","note":"blank=untested, NOT safe/additive"})
mat=pd.DataFrame(rows)
mat.to_parquet(ROOT/"results/tables/PAIRWISE_INTERACTION_MATRIX.parquet", index=False)
mat.to_csv(ROOT/"results/tables/PAIRWISE_INTERACTION_MATRIX.csv", index=False)
print(f"pairs: {len(mat)}")
print(mat["relation"].value_counts().to_string())
print("\nnon-trivial (redundant/shared) pairs:")
print(mat[mat.relation!="untested/unknown"][["intervention_a","intervention_b","relation"]].to_string(index=False))
