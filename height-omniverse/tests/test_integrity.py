"""Phase-A integrity QA (>=25 checks). Run: .venv/bin/python -m pytest tests/test_integrity.py -q"""
import duckdb, pathlib, re, pytest
ROOT = pathlib.Path(__file__).resolve().parents[1]
DB = ROOT / "database" / "height_omniverse.duckdb"

@pytest.fixture(scope="module")
def con():
    c = duckdb.connect(str(DB), read_only=True); yield c; c.close()
def one(con, q): return con.execute(q).fetchone()[0]

# --- provenance / source integrity ---
def test_migration_ref_tables_exist(con):
    have={r[0] for r in con.execute("SHOW TABLES").fetchall()}
    assert {"ref_species_fusion","ref_endpoint_rank"} <= have
def test_every_edge_has_source(con):
    assert one(con,"SELECT count(*) FROM evidence_edges WHERE source_id IS NULL OR source_id=''")==0
def test_every_calibration_row_has_source(con):
    assert one(con,"SELECT count(*) FROM surrogate_calibration WHERE source_id IS NULL OR source_id=''")==0
def test_every_study_has_source_location(con):
    assert one(con,"SELECT count(*) FROM studies WHERE source_location IS NULL OR source_location=''")==0
def test_every_intervention_has_provenance(con):
    assert one(con,"SELECT count(*) FROM interventions WHERE source_provenance IS NULL OR source_provenance=''")==0
def test_access_failures_have_route(con):
    assert one(con,"SELECT count(*) FROM access_failures WHERE request_route IS NULL OR request_route=''")==0

# --- endpoint hierarchy integrity ---
def test_studies_have_endpoint_rank(con):
    assert one(con,"SELECT count(*) FROM studies WHERE endpoint_rank IS NULL")==0
def test_final_length_flag_consistent(con):
    # anything flagged final-length must rank 1-2
    assert one(con,"SELECT count(*) FROM studies WHERE is_final_length=TRUE AND endpoint_rank>2")==0
def test_velocity_not_flagged_final(con):
    assert one(con,"SELECT count(*) FROM studies WHERE endpoint_type='velocity' AND is_final_length=TRUE")==0
def test_ref_endpoint_rank_bounds(con):
    assert one(con,"SELECT count(*) FROM ref_endpoint_rank WHERE endpoint_rank<1 OR endpoint_rank>10")==0

# --- species / fusion integrity ---
def test_studies_have_fusion_flag(con):
    assert one(con,"SELECT count(*) FROM studies WHERE fusion_capable_flag IS NULL OR fusion_capable_flag=''")==0
def test_rodent_not_fusion_capable(con):
    assert one(con,"SELECT count(*) FROM studies WHERE (lower(species) LIKE '%mouse%' OR lower(species) LIKE '%rat%') AND fusion_capable_flag='yes'")==0
def test_human_is_fusion_capable(con):
    assert one(con,"SELECT count(*) FROM studies WHERE lower(species) LIKE '%human%' AND fusion_capable_flag='no'")==0

# --- no recommendation / dose-as-recommendation leakage ---
def test_no_recommendation_columns(con):
    for t in ["interventions","studies","surrogate_calibration","evidence_edges"]:
        cols={r[0].lower() for r in con.execute(f"DESCRIBE {t}").fetchall()}
        assert not any("recommend" in c or c in ("best_stack","regimen","protocol") for c in cols), t
def test_no_first_person_recommendation(con):
    # WE never recommend; quoting a paper's clinical caution ("not routinely recommended") is allowed.
    bad=one(con,"SELECT count(*) FROM surrogate_calibration WHERE lower(notes) LIKE '%we recommend%' OR lower(notes) LIKE '%should take%' OR lower(notes) LIKE '%best stack%' OR lower(notes) LIKE '%optimal dose%'")
    assert bad==0

# --- simulation must not enter empirical tables ---
def test_no_simulation_in_studies(con):
    assert one(con,"SELECT count(*) FROM studies WHERE lower(design) LIKE '%simulation%' OR lower(intervention) LIKE '%simulation%'")==0
def test_simulation_edges_are_tier8(con):
    # any edge whose source is a model/simulation must be T8
    bad=one(con,"SELECT count(*) FROM evidence_edges WHERE (lower(source_id) LIKE '%simulation%' OR lower(source_id) LIKE '%model%' OR source_id='companion analysis') AND evidence_tier<>'T8'")
    assert bad==0

# --- dedup ---
def test_no_duplicate_gene_symbols(con):
    assert one(con,"SELECT count(*) FROM (SELECT symbol FROM genes GROUP BY symbol HAVING count(*)>1)")==0
def test_no_duplicate_intervention_ids(con):
    assert one(con,"SELECT count(*) FROM (SELECT intervention_id FROM interventions GROUP BY intervention_id HAVING count(*)>1)")==0
def test_no_duplicate_edge_ids(con):
    assert one(con,"SELECT count(*) FROM (SELECT edge_id FROM evidence_edges GROUP BY edge_id HAVING count(*)>1)")==0

# --- controlled vocabularies ---
def test_directness_tiers_valid(con):
    vals={r[0] for r in con.execute("SELECT DISTINCT directness_tier FROM studies").fetchall()}
    assert vals <= {"T1","T2","T3","T4","T5","T6","T7","T8"}
def test_edge_tiers_valid(con):
    vals={r[0] for r in con.execute("SELECT DISTINCT evidence_tier FROM evidence_edges").fetchall()}
    assert vals <= {"T1","T2","T3","T4","T5","T6","T7","T8"}
def test_direction_values_valid(con):
    vals={r[0] for r in con.execute("SELECT DISTINCT direction FROM evidence_edges").fetchall()}
    assert vals <= {"+","-","0","?","+/-"}

# --- gene universe / sanity ---
def test_gene_universe_size(con):
    assert one(con,"SELECT count(*) FROM genes")>40000
def test_pmid_format_where_present(con):
    rows=con.execute("SELECT source_id FROM evidence_edges WHERE source_id LIKE 'PMID%'").fetchall()
    assert all(re.search(r'PMID\s*\d+', r[0]) for r in rows)

# --- contradiction integrity ---
def test_contradiction_flag_populated(con):
    # letrozole PAH vs adult-height edges must be linked as contradictions
    assert one(con,"SELECT count(*) FROM evidence_edges WHERE has_contradiction=TRUE")>=2
def test_views_exist(con):
    have={r[0] for r in con.execute("SHOW TABLES").fetchall()}
    for v in ["v_adult_height_evidence","v_velocity_pos_final_null","v_brake_genes",
              "v_final_length_after_washout","v_mouse_or_condyle_only","v_contradictions"]:
        assert v in have, v
def test_velocity_pos_final_null_nonempty(con):
    assert one(con,"SELECT count(*) FROM v_velocity_pos_final_null")>=3
