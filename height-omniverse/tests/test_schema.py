"""Schema-integrity, provenance, and dedup QA tests. Run: .venv/bin/python -m pytest tests/ -q"""
import duckdb, pathlib, pytest
ROOT = pathlib.Path(__file__).resolve().parents[1]
DB = ROOT / "database" / "height_omniverse.duckdb"

@pytest.fixture(scope="module")
def con():
    c = duckdb.connect(str(DB), read_only=True); yield c; c.close()

def test_core_tables_exist(con):
    have = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
    for t in ["genes","interventions","studies","outcomes","evidence_edges",
              "evidence_snippets","access_failures","surrogate_calibration"]:
        assert t in have, f"missing table {t}"

def test_gene_universe_is_genome_wide(con):
    assert con.execute("SELECT count(*) FROM genes").fetchone()[0] > 40000

def test_every_edge_has_a_source(con):
    n = con.execute("SELECT count(*) FROM evidence_edges").fetchone()[0]
    if n == 0: pytest.skip("no edges loaded")
    bad = con.execute("SELECT count(*) FROM evidence_edges WHERE source_id IS NULL OR source_id=''").fetchone()[0]
    assert bad == 0, f"{bad} edges without source_id"

def test_every_study_has_directness_tier(con):
    n = con.execute("SELECT count(*) FROM studies").fetchone()[0]
    if n == 0: pytest.skip("no studies loaded")
    bad = con.execute("SELECT count(*) FROM studies WHERE directness_tier IS NULL OR directness_tier=''").fetchone()[0]
    assert bad == 0

def test_calibration_translation_values_controlled(con):
    n = con.execute("SELECT count(*) FROM surrogate_calibration").fetchone()[0]
    if n == 0: pytest.skip("no calibration rows")
    vals = {r[0] for r in con.execute("SELECT DISTINCT upper(surrogate_translated) FROM surrogate_calibration").fetchall()}
    allowed = {"YES","NO","UNKNOWN","PARTIAL","PARTIAL/MIXED","UNKNOWN-SUGGESTIVE",""}
    assert vals <= allowed, f"uncontrolled translation values: {vals - allowed}"

def test_no_duplicate_gene_symbols(con):
    dup = con.execute("SELECT count(*) FROM (SELECT symbol FROM genes GROUP BY symbol HAVING count(*)>1)").fetchone()[0]
    assert dup == 0, f"{dup} duplicate gene symbols"
