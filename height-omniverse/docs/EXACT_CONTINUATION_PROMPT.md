# Exact continuation prompt

Resume the Height Omniverse completion on branch `claude/height-omniverse-data-harvest` (PR #2).
State: Phases A, D, E, F, G, H, K, M, O complete or substantially done; B/C/N partial; I/J deferred.

## Resume commands
```bash
cd height-omniverse
bash run_all.sh                       # rebuild DB + models (needs GSE288529 h5ad for D/E; guarded)
.venv/bin/python -m pytest tests/ -q  # expect 34 passing
```

## Next executable tasks (in priority order), each ends in commit + PR update
1. **B (finish saturation):** run 2 more forward/backward citation rounds per category via
   `src/literature/*` (Europe PMC + NCBI elink); append to SCREENING_LOG/CITATION_FRONTIER until a full
   round yields no new included item; update PRISMA_COUNTS. Acceptance: 3 rounds/category logged.
2. **F2 (LINCS):** run `src/interventions/lincs_run.py` (to be written) against GSE92742/70138 subsets;
   build empirical growth-plate desired-state signature from GENOME_WIDE_ZONE_EXPRESSION.parquet;
   connectivity = hypothesis only. Mark partial if L1000 too large; run largest valid subset.
3. **D2 (omics from raw):** download GSE244884 / GSE275171 / GSE288028 processed matrices; harmonize
   into the zonal atlas; donor-level pseudobulk; keep species/age separate. Acceptance: each processed
   or reproducibly blocked.
4. **C1 (regulatory/full-text):** PMC XML for the 62 included PMIDs (table extraction); FDA/EMA PDFs via
   pdftotext where lawful; extract juvenile-tox/body-length/physeal endpoints. Locators required.
5. **I1 (agent-based lineage):** implement `src/modeling/abm_lineage.py` (stochastic R0→R1→P→PH→H with
   neutral clone competition); calibrate to lineage-tracing clone-size distributions.
6. **J1 (FE):** if FEniCSx installable, implement zonal mechanobiology; else keep spatial RD surrogate.
7. **N1 (cards):** enrich `src/literature/make_cards.py` to embed zonal topology + compounds +
   literature + contradictions; empty fields say "no qualifying evidence located after N rounds".

## Invariants to keep
Endpoint hierarchy; no recommendation/dose fields; simulation labeled simulation; cells ≠ biological
replicates; every load-bearing value has a source_id; update COMPLETION_BACKLOG.csv + COMPLETION_STATE.md
after each phase; commit + push + update PR #2.
