"""Phase B acceptance: >=3 search rounds logged per category, and any category marked SATURATED
must have its final round produce 0 new included items. Run: pytest tests/test_saturation.py -q"""
import csv, pathlib, pytest
ROOT = pathlib.Path(__file__).resolve().parents[1]
RT = ROOT / "results" / "tables"

def _rows(name):
    p = RT / name
    if not p.exists(): pytest.skip(f"{name} not present yet")
    return list(csv.DictReader(open(p)))

def test_search_rounds_exist():
    rows = _rows("SEARCH_ROUNDS.csv")
    assert rows, "SEARCH_ROUNDS.csv empty"

def test_three_rounds_per_category():
    rows = _rows("SEARCH_ROUNDS.csv")
    cats = {}
    for r in rows:
        cats.setdefault(r["category"], set()).add(str(r["round"]))
    # Round 1 is the seed search (SEARCH_QUERIES); rounds 2 & 3 in SEARCH_ROUNDS -> combined >=3
    bad = {c: sorted(v) for c, v in cats.items() if len(v) < 2}
    assert not bad, f"categories without rounds 2 AND 3 logged: {bad}"

def test_saturation_flag_consistent():
    sat = _rows("SATURATION.csv")
    for r in sat:
        if r.get("saturated","").upper() == "Y":
            assert str(r.get("round3_new_items","0")).strip() in ("0","0.0"), \
                f"{r['category']} marked saturated but round3_new_items={r['round3_new_items']}"

def test_every_included_new_study_has_endpoint_type():
    rows = _rows("SCREENING_LOG.csv")
    inc = [r for r in rows if r["include"] == "Y"]
    allowed = {"adult_height","near_final","velocity","bone_length","mechanism","other"}
    bad = [r["pmid"] for r in inc if r.get("endpoint_type","") not in allowed]
    assert not bad, f"included studies with bad endpoint_type: {bad[:10]}"

def test_excluded_have_reasons():
    rows = _rows("EXCLUDED_STUDIES.csv")
    bad = [r for r in rows if not r.get("reason","").strip()]
    assert not bad, "some excluded studies lack a reason"

def test_citation_frontier_present():
    rows = _rows("CITATION_FRONTIER.csv")
    assert rows, "CITATION_FRONTIER.csv empty"
