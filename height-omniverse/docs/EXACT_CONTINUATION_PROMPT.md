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

## UPDATE (post rounds 2-3): Phase B round 4 for 11 unsaturated categories
Saturated (stop): Q1, Q3, Q13, Q14. Run round 4 (forward+backward elink) for the rest, seeded on the
round-3 newly-included PMIDs; a category is done when a round yields 0 new. Round-4 seeds:
- Q2: 27710241,40329244 | Q4: 30335613,34341520 | Q5: 30133462,39571550 | Q6: 15784709,20061421
- Q7: 25824333,29073591 | Q8: 37080994,39337434 | Q9: 19730156,16614378 | Q10: 32282051,30235256
- Q11: 39423254,39991457 | Q12: 38027309,22489887 | Q15: 35330881,36387899
Dedup against results/tables/SCREENING_LOG.csv (212 PMIDs). Note: Q8/Q10/Q15 are deep mechanism/genetics
citation networks and may need several more rounds; log honestly, never fake saturation.
Append to SEARCH_ROUNDS/CITATION_FRONTIER/SCREENING_LOG/EXCLUDED_STUDIES/SATURATION, rebuild, re-run
pytest (40+ pass), re-run make_cards, update PRISMA + backlog, commit, push, update PR #2.
