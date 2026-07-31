#!/usr/bin/env python
"""Build the Height Omniverse DuckDB warehouse: normalized schema + genome-wide gene universe.
Idempotent: DROP+CREATE tables, reload from data/raw + curated CSVs. Never overwrites raw files.
Every load-bearing row must trace to a source_id."""
import duckdb, pandas as pd, pathlib, hashlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DB = ROOT / "database" / "height_omniverse.duckdb"
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
CUR = ROOT / "data" / "curated"
PROC.mkdir(exist_ok=True, parents=True)

con = duckdb.connect(str(DB))

# ---------------- schema ----------------
con.execute("""
CREATE OR REPLACE TABLE genes (
  ensembl_id VARCHAR, hgnc_id VARCHAR, symbol VARCHAR, name VARCHAR,
  locus_group VARCHAR, locus_type VARCHAR, biotype VARCHAR, chromosome VARCHAR,
  aliases VARCHAR, prev_symbols VARCHAR, gene_group VARCHAR, entrez_id VARCHAR, status VARCHAR,
  height_gwas_evidence VARCHAR, monogenic_stature_evidence VARCHAR,
  growthplate_zone_evidence VARCHAR, direction_of_effect VARCHAR, source_provenance VARCHAR
);
CREATE OR REPLACE TABLE interventions (
  intervention_id VARCHAR, name VARCHAR, synonyms VARCHAR, modality VARCHAR, target VARCHAR,
  action VARCHAR, approval_status VARCHAR, route VARCHAR, pk_summary VARCHAR,
  evidence_context VARCHAR, technology_readiness VARCHAR, source_provenance VARCHAR
);
CREATE OR REPLACE TABLE studies (
  study_id VARCHAR, citation VARCHAR, year INTEGER, design VARCHAR, species VARCHAR, sex VARCHAR,
  age_stage VARCHAR, disease VARCHAR, tissue_site VARCHAR, intervention VARCHAR, comparator VARCHAR,
  n_true_bioreplicates VARCHAR, exposure VARCHAR, duration VARCHAR, washout VARCHAR, endpoint VARCHAR,
  maturity_followup VARCHAR, effect_size VARCHAR, uncertainty VARCHAR, direction VARCHAR,
  geometry_mineralization VARCHAR, directness_tier VARCHAR, risk_of_bias VARCHAR, source_location VARCHAR, notes VARCHAR
);
CREATE OR REPLACE TABLE outcomes (
  outcome_id VARCHAR, study_id VARCHAR, outcome_type VARCHAR, raw_value VARCHAR, standardized_value VARCHAR,
  unit VARCHAR, endpoint_rank INTEGER, is_final_length BOOLEAN, notes VARCHAR
);
CREATE OR REPLACE TABLE evidence_edges (
  edge_id VARCHAR, subject VARCHAR, predicate VARCHAR, object VARCHAR, direction VARCHAR, context VARCHAR,
  effect_size VARCHAR, source_id VARCHAR, evidence_tier VARCHAR, directness VARCHAR, confidence VARCHAR,
  conflicting_edges VARCHAR
);
CREATE OR REPLACE TABLE evidence_snippets (
  snippet_id VARCHAR, source_id VARCHAR, excerpt VARCHAR, paraphrase VARCHAR, location VARCHAR,
  extractor_confidence VARCHAR, human_reviewed BOOLEAN, links_to VARCHAR
);
CREATE OR REPLACE TABLE access_failures (
  source VARCHAR, reason VARCHAR, attempted_method VARCHAR, metadata_found VARCHAR,
  request_route VARCHAR, importance VARCHAR
);
CREATE OR REPLACE TABLE literature (
  study_id VARCHAR, pmid VARCHAR, title VARCHAR, year VARCHAR, journal VARCHAR,
  query_id VARCHAR, endpoint_type VARCHAR, source VARCHAR
);
CREATE OR REPLACE TABLE compound_target (
  target_gene VARCHAR, molecule_chembl_id VARCHAR, pref_name VARCHAR, max_phase VARCHAR,
  standard_type VARCHAR, standard_value_nM VARCHAR, source VARCHAR
);
CREATE OR REPLACE TABLE surrogate_calibration (
  intervention VARCHAR, population VARCHAR, design VARCHAR, n VARCHAR, surrogate_effect VARCHAR,
  adult_height_effect VARCHAR, surrogate_translated VARCHAR, direction VARCHAR, source_id VARCHAR, notes VARCHAR
);
""")

# ---------------- genome-wide gene universe (HGNC complete set) ----------------
hgnc = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False).fillna("")
g = pd.DataFrame({
    "ensembl_id": hgnc["ensembl_gene_id"], "hgnc_id": hgnc["hgnc_id"], "symbol": hgnc["symbol"],
    "name": hgnc["name"], "locus_group": hgnc["locus_group"], "locus_type": hgnc["locus_type"],
    "biotype": hgnc["locus_type"], "chromosome": hgnc["location"],
    "aliases": hgnc["alias_symbol"], "prev_symbols": hgnc["prev_symbol"], "gene_group": hgnc["gene_group"],
    "entrez_id": hgnc["entrez_id"], "status": hgnc["status"],
    "height_gwas_evidence": "", "monogenic_stature_evidence": "", "growthplate_zone_evidence": "",
    "direction_of_effect": "", "source_provenance": "HGNC complete set (retrieved 2026-07-31)",
})
g = g[g["status"] == "Approved"].reset_index(drop=True)
con.execute("DELETE FROM genes")
con.register("g_df", g)
con.execute("INSERT INTO genes SELECT * FROM g_df")
g.to_parquet(PROC / "ALL_GENES.parquet", index=False)

# ---------------- load curated CSVs if present (created by ingest_curated.py) ----------------
def load_csv(table, fname):
    p = CUR / fname
    if p.exists():
        df = pd.read_csv(p, dtype=str).fillna("")
        con.register("df_tmp", df)
        cols = [r[0] for r in con.execute(f"DESCRIBE {table}").fetchall()]
        common = [c for c in cols if c in df.columns]
        con.execute(f"INSERT INTO {table} ({','.join(common)}) SELECT {','.join(common)} FROM df_tmp")
        con.unregister("df_tmp")
        return len(df)
    return 0

loaded = {t: load_csv(t, f) for t, f in [
    ("interventions","interventions.csv"), ("studies","studies.csv"), ("outcomes","outcomes.csv"),
    ("evidence_edges","evidence_edges.csv"), ("surrogate_calibration","surrogate_calibration.csv"),
    ("access_failures","access_failures.csv"), ("evidence_snippets","evidence_snippets.csv"),
    ("literature","literature.csv"), ("compound_target","compound_target.csv")]}

# annotate gene flags from curated genetics table if present
gp = CUR / "genetic_direction.csv"
if gp.exists():
    gd = pd.read_csv(gp, dtype=str).fillna("")
    for _, r in gd.iterrows():
        con.execute("UPDATE genes SET monogenic_stature_evidence=?, direction_of_effect=?, source_provenance=source_provenance||'; '||? WHERE symbol=?",
                    [r.get("key_phenotype",""), r.get("height_direction",""), r.get("source_id",""), r.get("gene","")])

# ---------------- apply migration_002 (additive schema hardening) ----------------
for mfn in ["migration_002.sql", "migration_003_views.sql"]:
    mig = ROOT / "database" / "schemas" / mfn
    if mig.exists():
        con.execute(mig.read_text()); print(f"applied {mfn}")

n_genes = con.execute("SELECT count(*) FROM genes").fetchone()[0]
print(f"genes (Approved): {n_genes}")
print("curated loaded:", {k:v for k,v in loaded.items() if v})
for t in ["genes","interventions","studies","outcomes","evidence_edges","surrogate_calibration","access_failures"]:
    print(f"  {t:22s}: {con.execute(f'SELECT count(*) FROM {t}').fetchone()[0]} rows")
con.close()
print("DB built:", DB)
