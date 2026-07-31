# Completion state

- **Branch:** claude/height-omniverse-data-harvest
- **Completed phase:** A (schema hardening + QA)
- **Next phase:** D (genome-wide zonal atlas) + E (GWAS×zone enrichment); B/F harvest agents running.
- **Resume command:** `cd height-omniverse && bash run_all.sh && .venv/bin/python -m pytest tests/ -q`
- **Storage:** DB rebuildable (gitignored); HGNC raw ~15 MB; GSE288529 processed reused from ../abaloparatide-growth.
- **Datasets downloaded:** HGNC complete set (checksummed). GSE288529 processed h5ad reused.
- **Unsuccessful URLs/APIs:** none this phase.
- **Pending acceptance tests:** B (3 saturation rounds), D (genome-wide parquet), E (enrichment), F (ChEMBL), H (spatial model mass-conservation).
- **Tests:** 34 passing (6 original + 28 integrity).
