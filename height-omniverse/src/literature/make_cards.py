#!/usr/bin/env python
"""Phase N: generate machine (JSON) + human (MD) evidence cards for every intervention and every
gene with recorded evidence. Cards embed genetics, zonal topology, compounds, literature, interactions,
and contradictions. Empty fields say 'no qualifying evidence located after N completed search rounds'.
NO recommendation field (by design). Every load-bearing line points to a source."""
import duckdb, pathlib, json, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "evidence_cards"; OUT.mkdir(parents=True, exist_ok=True)
TB = ROOT / "results" / "tables"
con = duckdb.connect(str(ROOT / "database" / "height_omniverse.duckdb"))
def rows(q, p=None): return con.execute(q, p or []).fetchdf().to_dict("records")

# completed search rounds (round 1 seed + any rounds in SEARCH_ROUNDS)
NROUNDS = 1
srp = TB / "SEARCH_ROUNDS.csv"
if srp.exists():
    import csv
    rs = {r["round"] for r in csv.DictReader(open(srp))}
    NROUNDS = len(rs)
EMPTY = f"no qualifying evidence located after {NROUNDS} completed search round(s)"

spec = pd.read_parquet(TB/"GENOME_WIDE_STATE_SPECIFICITY.parquet") if (TB/"GENOME_WIDE_STATE_SPECIFICITY.parquet").exists() else None
dual = pd.read_csv(TB/"DUAL_EVIDENCE_CANDIDATES.csv") if (TB/"DUAL_EVIDENCE_CANDIDATES.csv").exists() else None
def mouse(s): return s[:1].upper()+s[1:].lower()

cards = 0
# ---- intervention cards ----
for iv in rows("SELECT * FROM interventions"):
    name = iv["name"]; key = name.split()[0]
    cal = rows("SELECT * FROM surrogate_calibration WHERE intervention ILIKE ?", [f"%{key}%"])
    edg = rows("SELECT * FROM evidence_edges WHERE subject ILIKE ?", [f"%{key}%"])
    inter = rows("SELECT intervention_a,intervention_b,relation FROM read_csv_auto('results/tables/PAIRWISE_INTERACTION_MATRIX.csv') WHERE (intervention_a=? OR intervention_b=?) AND relation<>'untested/unknown'", [name,name])
    card = {"type":"intervention","name":name,"modality":iv.get("modality",""),"target":iv.get("target",""),
            "action":iv.get("action",""),"approval_status":iv.get("approval_status",""),"route":iv.get("route",""),
            "pk":iv.get("pk_summary",""),"surrogate_vs_final_height":cal or EMPTY,"evidence_edges":edg or EMPTY,
            "interactions":inter or EMPTY,"safety_note":"efficacy/safety stored separately; see EVIDENCE_LEDGER (companion)",
            "recommendation":"NONE BY DESIGN","source_provenance":iv.get("source_provenance","")}
    slug="".join(c if c.isalnum() else "_" for c in name)[:60]
    (OUT/f"intervention__{slug}.json").write_text(json.dumps(card,indent=2,default=str))
    md=[f"# {name} — intervention evidence card","",f"- **Modality:** {iv.get('modality','')}  ",
        f"- **Target/action:** {iv.get('target','')} / {iv.get('action','')}  ",
        f"- **Approval:** {iv.get('approval_status','')}  · **Route:** {iv.get('route','')}  · **PK:** {iv.get('pk_summary','')}","",
        "## Surrogate vs ATTAINED final height"]
    md += [f"- {c['population']} ({c['design']}): surrogate={c['surrogate_effect']}; adult height={c['adult_height_effect']}; **translated={c['surrogate_translated']}** [{c['source_id']}]" for c in cal] or [f"- {EMPTY}"]
    md += ["","## Interactions (pathway convergence; blank=untested, NOT safe/additive)"]
    md += [f"- {i['intervention_a']} × {i['intervention_b']}: {i['relation']}" for i in inter] or [f"- {EMPTY}"]
    md += ["", f"_Provenance: {iv.get('source_provenance','')}_", "", "> No recommendation field by design."]
    (OUT/f"intervention__{slug}.md").write_text("\n".join(md)); cards+=2

# ---- gene cards (genes with recorded evidence) ----
for g in rows("SELECT * FROM genes WHERE monogenic_stature_evidence <> '' OR direction_of_effect <> ''"):
    sym=g["symbol"]; m=mouse(sym)
    topo = None
    if spec is not None:
        sr = spec[spec["gene"]==m]
        if len(sr): topo={"top_zone":sr["max_state"].iloc[0],"specificity":round(float(sr["specificity_score"].iloc[0]),2),
                          "gp_expression":round(float(sr[[f"mean_{s}" for s in ['Resting','Proliferative','Prehypertrophic','Hypertrophic']]].sum(1).iloc[0]),2)}
    comp = rows("SELECT molecule_chembl_id,pref_name,max_phase,standard_type,standard_value_nM FROM compound_target WHERE upper(target_gene)=?", [sym.upper()])
    lit = rows("SELECT pmid,title,year FROM literature WHERE upper(title) LIKE ?", [f"%{sym.upper()}%"])
    de = dual[dual["gene"]==sym].to_dict("records")[0] if (dual is not None and (dual["gene"]==sym).any()) else None
    card={"type":"gene","symbol":sym,"name":g.get("name",""),"ensembl":g.get("ensembl_id",""),
          "genetics_direction":g.get("direction_of_effect",""),"monogenic_phenotype":g.get("monogenic_stature_evidence","") or EMPTY,
          "zonal_topology":topo or f"{EMPTY} (or absent in mouse atlas)","compounds":comp or EMPTY,
          "literature_mentions":lit or EMPTY,"dual_evidence":de or EMPTY,
          "contradictions":"see evidence_edges conflicting_edges","recommendation":"NONE BY DESIGN",
          "source_provenance":g.get("source_provenance","")}
    (OUT/f"gene__{sym}.json").write_text(json.dumps(card,indent=2,default=str))
    comp_str = ", ".join(f"{c['pref_name'] or c['molecule_chembl_id']}({c['standard_type']} {c['standard_value_nM']}nM)" for c in comp[:6]) if comp else EMPTY
    de_str = (f"{de['n_evidence_modes']}" + (" (high-priority)" if de.get('high_priority_research') else "")) if de else EMPTY
    md=[f"# {sym} — gene evidence card  ({g.get('name','')})","",
        f"- **Genetics direction:** {g.get('direction_of_effect','')} — {g.get('monogenic_stature_evidence','') or EMPTY}  ",
        f"- **Zonal topology:** {topo or EMPTY}  ",
        f"- **Compounds (ChEMBL):** {comp_str}",
        f"- **Dual-evidence modes:** {de_str}","",
        f"_Provenance: {g.get('source_provenance','')}_","", "> No recommendation field by design."]
    (OUT/f"gene__{sym}.md").write_text("\n".join(md)); cards+=2
con.close()
print(f"wrote {cards} evidence-card files (rounds completed: {NROUNDS}); empties labeled '{EMPTY}'")
