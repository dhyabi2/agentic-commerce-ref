# Limitations — what is scaffolded vs run

This is a **real reproducible subset**, not the full multi-month spec. Honest status:

## Built and run
DuckDB warehouse + 45,030-gene universe; ontologies; surrogate→final-height calibration (17 rows);
intervention atlas (14, verified PMIDs); genetics direction-of-effect (24 rows, verified PMIDs);
mechanical + nutrient/botanical ledgers; omics-availability catalogue (11 datasets); knowledge-graph
GraphML; 49 evidence cards; seed digital-twin (reserve-limited growth-plate ODE + PK/PD receptor
reserve, reused from the companion abaloparatide analysis); QA tests (6 passing); dashboard.

## Scaffolded / documented but NOT run to completion (do not treat as complete)
- **Genome-wide omics reprocessing from raw** — only the prior GSE288529 single-library audit is
  reused; the other 10+ catalogued datasets are availability-logged, not reprocessed.
- **LINCS/CMap, ChEMBL/PubChem/BindingDB** compound & perturbation layers — not ingested.
- **GWAS fine-mapping / MR / LDSC-SEG cell-type enrichment** — GWAS captured at anchor level only;
  the growth-plate cell-type enrichment narrative is explicitly UNVERIFIED for Yengo 2022.
- **FEniCS finite-element mechanical twin, spatial reaction-diffusion, agent-based lineage models** —
  the ODE + PK/PD twin runs; the other model classes are specified, not implemented.
- **Neo4j server / vector full-text index** — GraphML export provided; no running Neo4j; the DuckDB +
  Parquet layer is the query surface.
- **Recursive citation-closure to formal saturation** — one seed round + targeted expansion; saturation
  criterion not met (see COVERAGE_REPORT).
- **Patents / theses / conference deep-mining** — partial.

## Epistemic caveats (enforced)
Simulation is labeled simulation; cells are never biological replicates; mRNA≠protein≠activity;
surrogate≠final height; species-fusion flags required; efficacy and safety kept separate; no dosing,
timing, sourcing, or "best stack" output anywhere.
