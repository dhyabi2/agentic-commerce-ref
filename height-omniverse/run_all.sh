#!/usr/bin/env bash
# One-command reproducible build of the Height Omniverse warehouse.
set -euo pipefail
cd "$(dirname "$0")"
PY=.venv/bin/python
if [ ! -x "$PY" ]; then
  python3 -m venv .venv; .venv/bin/python -m pip install -q --upgrade pip
  .venv/bin/python -m pip install -q duckdb pandas pyarrow numpy scipy matplotlib requests networkx anndata h5py pytest
fi
if [ ! -s data/raw/hgnc_complete_set.txt ]; then
  echo "== fetch HGNC gene universe =="
  curl -sSL --retry 4 -o data/raw/hgnc_complete_set.txt \
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
fi
echo "== [1] curated verified evidence + compounds + literature =="; $PY src/ingestion/build_curated.py
# --- omics/discovery steps that need the external GSE288529 processed h5ad (Phase D/E) ---
H5=/home/user/agentic-commerce/growth_plate_pth1r/data/processed.h5ad
if [ -f "$H5" ]; then
  echo "== [2] Phase D genome-wide zonal atlas =="; $PY src/omics/genome_wide_zone_expression.py
  echo "== [3] Phase E GWAS x zone enrichment =="; $PY src/genetics/gwas_zone_enrichment.py
else
  echo "== [2-3] SKIPPED (GSE288529 processed h5ad absent; see docs/BLOCKED_DATA.md) =="
fi
echo "== [4] Phase K interaction matrix =="; $PY src/interventions/build_interactions.py
echo "== [5] Phase G dual-evidence discovery =="; $PY src/genetics/dual_evidence.py || true
echo "== [6] build DuckDB (schema + migrations + views + all layers) =="; $PY src/ingestion/build_db.py
echo "== [7] knowledge graph =="; $PY src/qa/export_graph.py
echo "== [8] evidence cards =="; $PY src/literature/make_cards.py
echo "== [9] Phase H spatial reaction-diffusion (SIM) =="; $PY src/mechanics/spatial_rd.py
echo "== [10] digital twin (PK/PD + growth-plate ODE) (SIM) =="; $PY src/pkpd/pkpd_receptor_reserve.py; $PY src/modeling/twin_intervention_sweep.py
echo "== [11] QA tests =="; $PY -m pytest tests/ -q
echo "DONE. Query DB: .venv/bin/python -c \"import duckdb;print(duckdb.connect('database/height_omniverse.duckdb').sql('SHOW TABLES'))\""
