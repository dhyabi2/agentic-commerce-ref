#!/usr/bin/env python
"""Generate machine-readable (JSON) + human-readable (MD) evidence cards for every intervention
and every gene with recorded evidence. No recommendation field (by design)."""
import duckdb, pathlib, json
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "evidence_cards"; OUT.mkdir(parents=True, exist_ok=True)
con = duckdb.connect(str(ROOT / "database" / "height_omniverse.duckdb"))

def rows(q, p=None): return con.execute(q, p or []).fetchdf().to_dict("records")

cards = 0
# intervention cards
for iv in rows("SELECT * FROM interventions"):
    name = iv["name"]
    cal = rows("SELECT * FROM surrogate_calibration WHERE intervention ILIKE ?", [f"%{name.split()[0]}%"])
    edg = rows("SELECT * FROM evidence_edges WHERE subject ILIKE ?", [f"%{name.split()[0]}%"])
    card = {"type": "intervention", "name": name, "aliases": iv.get("synonyms",""),
            "modality": iv.get("modality",""), "target": iv.get("target",""),
            "action": iv.get("action",""), "approval_status": iv.get("approval_status",""),
            "route": iv.get("route",""), "pk": iv.get("pk_summary",""),
            "surrogate_vs_final_height": cal, "evidence_edges": edg,
            "source_provenance": iv.get("source_provenance","")}
    slug = "".join(c if c.isalnum() else "_" for c in name)[:60]
    (OUT / f"intervention__{slug}.json").write_text(json.dumps(card, indent=2, default=str))
    md = [f"# {name}  (intervention evidence card)", "",
          f"- **Modality:** {iv.get('modality','')}", f"- **Target:** {iv.get('target','')}",
          f"- **Action:** {iv.get('action','')}", f"- **Approval:** {iv.get('approval_status','')}",
          f"- **Route:** {iv.get('route','')}", f"- **PK:** {iv.get('pk_summary','')}", "",
          "## Surrogate vs final height"]
    for c in cal:
        md.append(f"- **{c['population']}** ({c['design']}): surrogate={c['surrogate_effect']}; "
                  f"adult height={c['adult_height_effect']}; **translated={c['surrogate_translated']}** "
                  f"[{c['source_id']}]")
    if not cal: md.append("- (no surrogate/final-height row recorded)")
    md += ["", f"_Source provenance: {iv.get('source_provenance','')}_",
           "", "> No recommendation field by design (research warehouse)."]
    (OUT / f"intervention__{slug}.md").write_text("\n".join(md)); cards += 2

# gene cards (only genes with recorded evidence)
for g in rows("SELECT * FROM genes WHERE monogenic_stature_evidence <> '' OR direction_of_effect <> ''"):
    sym = g["symbol"]
    card = {"type": "gene", "symbol": sym, "name": g.get("name",""), "ensembl": g.get("ensembl_id",""),
            "monogenic_stature_evidence": g.get("monogenic_stature_evidence",""),
            "direction_of_effect": g.get("direction_of_effect",""),
            "growthplate_zone_evidence": g.get("growthplate_zone_evidence",""),
            "source_provenance": g.get("source_provenance","")}
    (OUT / f"gene__{sym}.json").write_text(json.dumps(card, indent=2, default=str)); cards += 1

con.close()
print(f"wrote {cards} evidence-card files to {OUT}")
