#!/usr/bin/env bash
# One-command reproducible pipeline for the HDAC4 growth plate analysis.
#
#   ./run_all.sh            full run (downloads ~760 MB from GEO on first use)
#   PY=/path/to/python ./run_all.sh    use a specific interpreter
#
# Environment is pinned in requirements.txt. Seeds are fixed in src/hdac4_config.py.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PY:-python3}"
VENV="${VENV:-.venv}"

if [ ! -x "$VENV/bin/python" ]; then
  echo "==> creating virtualenv in $VENV"
  "$PY" -m venv "$VENV"
  "$VENV/bin/pip" install --quiet --upgrade pip setuptools wheel
  "$VENV/bin/pip" install --quiet -r requirements.txt
fi
PYBIN="$VENV/bin/python"

mkdir -p data/raw data/processed results/tables results/figures docs

fetch () { # url dest
  local url="$1" dest="$2"
  [ -s "$dest" ] && { echo "    have $(basename "$dest")"; return 0; }
  for i in 1 2 3 4 5; do
    if curl -sSL -m 2400 "$url" -o "$dest"; then
      case "$dest" in
        *.tar) tar -tf "$dest" >/dev/null 2>&1 && return 0 ;;
        *.gz)  gzip -t "$dest" >/dev/null 2>&1 && return 0 ;;
        *)     return 0 ;;
      esac
    fi
    echo "    retry $i for $(basename "$dest")"; sleep $((2 ** i))
  done
  echo "FAILED to fetch $url" >&2; return 1
}

echo "==> 0. fetching data from GEO"
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288028&format=file" data/raw/GSE288028_RAW.tar
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE288529&format=file" data/raw/GSE288529_RAW.tar
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16981/matrix/GSE16981_series_matrix.txt.gz" data/raw/GSE16981_series_matrix.txt.gz
fetch "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1355/annot/GPL1355.annot.gz" data/raw/GPL1355.annot.gz
mkdir -p data/raw/GSE288028 data/raw/GSE288529
tar -xf data/raw/GSE288028_RAW.tar -C data/raw/GSE288028
tar -xf data/raw/GSE288529_RAW.tar -C data/raw/GSE288529
( cd data/raw && sha256sum *.tar *.gz ) > results/tables/raw_checksums.txt

echo "==> 1. GEO search (query log + hit counts)"
"$PYBIN" src/01_geo_search.py

echo "==> 2. load, QC, cluster, assign zones"
"$PYBIN" src/02_process.py

echo "==> 3. HDAC4 activity scores (decoupler ULM + scanpy score_genes)"
"$PYBIN" src/03_score.py

echo "==> 4. analyses A-D and negative controls"
"$PYBIN" src/04_analyses.py

echo "==> 5. supplementary rat temporal analysis (GSE16981)"
"$PYBIN" src/06_rat_temporal.py

echo "==> 6. figures"
"$PYBIN" src/05_figures.py

echo "==> 7. environment record"
"$PYBIN" -c "import json,sys,platform;from importlib.metadata import version,PackageNotFoundError
pkgs=['scanpy','anndata','decoupler','scikit-learn','pandas','numpy','statsmodels','scipy','matplotlib','h5py','leidenalg','igraph','numba']
d={'python':sys.version,'platform':platform.platform()}
for p in pkgs:
    try: d[p]=version(p)
    except PackageNotFoundError: d[p]='NOT INSTALLED'
open('results/tables/environment.json','w').write(json.dumps(d,indent=2))
print(json.dumps(d,indent=2))"

echo "==> done. See RESULTS.md, results/tables/ and results/figures/"
