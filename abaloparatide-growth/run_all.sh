#!/usr/bin/env bash
# One-command reproducible pipeline. Requires: python3 venv at .venv with numpy/scipy/pandas/matplotlib
# (see environment.txt). No network access needed — the receptor-topology numbers are the audited
# GSE288529 pseudobulk from the prior analysis (single library; descriptive only).
set -euo pipefail
cd "$(dirname "$0")"
PY=./.venv/bin/python
if [ ! -x "$PY" ]; then
  echo "creating venv..."; python3 -m venv .venv
  ./.venv/bin/python -m pip install -q --upgrade pip
  ./.venv/bin/python -m pip install -q numpy scipy pandas matplotlib
fi
echo "== [1/6] evidence ledger =="        ; $PY src/build_evidence_ledger.py
echo "== [2/6] claims audit =="           ; $PY src/build_claims_audit.py
echo "== [3/6] PK/PD receptor reserve ==" ; $PY src/pkpd_receptor_reserve.py
echo "== [4/6] growth-plate validation ==" ; $PY src/growth_plate_model.py
echo "== [5/6] growth-plate MC synthesis ==" ; $PY src/run_growth_synthesis.py
echo "== [6/6] evidence figures =="       ; $PY src/make_evidence_figs.py
echo "DONE. See REPORT.md, EXECUTIVE_SUMMARY.md, results/{tables,figures,models,sensitivity}/"
