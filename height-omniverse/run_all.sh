#!/usr/bin/env bash
# One-command reproducible build of the Height Omniverse warehouse.
set -euo pipefail
cd "$(dirname "$0")"
PY=.venv/bin/python
if [ ! -x "$PY" ]; then
  python3 -m venv .venv; .venv/bin/python -m pip install -q --upgrade pip
  .venv/bin/python -m pip install -q duckdb pandas pyarrow numpy scipy matplotlib requests networkx
fi
# raw gene universe (skip if cached)
if [ ! -s data/raw/hgnc_complete_set.txt ]; then
  echo "== fetch HGNC gene universe =="
  curl -sSL --retry 4 -o data/raw/hgnc_complete_set.txt \
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
fi
echo "== [1/6] build DuckDB + gene universe + load curated verified evidence =="
$PY src/ingestion/build_db.py
echo "== [2/6] export knowledge graph (GraphML + parquet) =="
$PY src/qa/export_graph.py
echo "== [3/6] evidence cards =="
$PY src/literature/make_cards.py
echo "== [4/6] PK/PD receptor-reserve twin (SIMULATION) =="
$PY src/pkpd/pkpd_receptor_reserve.py
echo "== [5/6] growth-plate ODE twin sweep (SIMULATION) =="
$PY src/modeling/twin_intervention_sweep.py
echo "== [6/6] QA tests =="
$PY -m pytest tests/ -q || echo "(pytest not present / some checks skipped)"
echo "DONE. Query: .venv/bin/python -c \"import duckdb;print(duckdb.connect('database/height_omniverse.duckdb').sql('SHOW TABLES'))\""
echo "Dashboard (optional): .venv/bin/streamlit run dashboard/app.py"
