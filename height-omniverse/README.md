# Height Omniverse — evidence warehouse & model-ready growth atlas

> **This repository is an evidence warehouse and simulation framework. It does not recommend a
> height intervention or personal protocol.**

A reproducible, queryable research corpus for one scientific question:

> *What biological, pharmacological, nutritional, mechanical, genetic, cellular, surgical, and
> tissue-engineering interventions can change stable final adult stature or mature segmental
> skeletal length, through what mechanisms, in which developmental states, and with what evidence?*

It is **not** a persuasive essay, **not** a therapeutic conclusion, and contains **no** dosing,
timing, sourcing, or "best stack" output. Experimental doses/PK are stored only as scientific
metadata (see `docs/METHODS.md` §research-only).

## Scope actually built in this session (honest)

This is a **real, reproducible subset** of the full spec, not the entire multi-month program. What
exists here is genuine and verified; what is not built is logged explicitly, not faked.

**Built & queryable:**
- DuckDB warehouse (`database/height_omniverse.duckdb`) with normalized schema + a genome-wide gene
  universe (**45,030 HGNC genes**, `data/processed/ALL_GENES.parquet`).
- Controlled-vocabulary **ontologies** (endpoint directness, species-fusion registry, developmental
  & cell states, 8-tier evidence hierarchy) — so evidence cannot be mixed accidentally.
- **Surrogate → final-height calibration** table (the crown jewel: velocity/predicted-height vs
  *attained* adult height, with verified PMIDs).
- **Intervention atlas** across all major modality classes, tiered, final-height separated from
  surrogates; **genetic direction-of-effect** table; **mechanical** & **nutrient/botanical** ledgers.
- Seed **digital twin** (reserve-limited growth-plate ODE + PK/PD receptor-reserve), reused from the
  companion abaloparatide analysis and clearly labeled **simulation, not evidence**.
- **Knowledge-graph** export (GraphML), evidence cards, coverage/blocked-data/VOI docs, example
  DuckDB queries, and a minimal dashboard.

**Explicitly NOT built (logged in `docs/BLOCKED_DATA.md`, `docs/LIMITATIONS.md`, `docs/COVERAGE_REPORT.md`):**
genome-wide reprocessing of every omics dataset from raw; full LINCS/CMap ingestion; FEniCS finite-
element mechanical twin; Neo4j server; recursive citation-closure to formal saturation; restricted
individual-level data (UKB/dbGaP/EGA). These are scaffolded and/or documented as data requests, not
claimed as complete.

## Reproduce
```bash
bash run_all.sh     # builds DB, loads curated verified evidence, runs models, exports graph + cards
```
Pins in `environment/python.lock`; seed 20260731; raw checksums in `docs/SOURCE_MANIFEST.csv`.

## Interrogate
See `docs/HOW_TO_INTERROGATE_THE_CORPUS.md` for example DuckDB queries (e.g. *all interventions with
human adult-height evidence*, *all candidates with positive velocity but null final-height*). Start
the dashboard with `./.venv/bin/streamlit run dashboard/app.py` (optional; requires streamlit).

## The one finding the corpus keeps surfacing (as data, not advice)
Across interventions with *both* surrogate and final-height data, **surrogate gains frequently do not
become adult-height gains** — GH banks height; aromatase inhibitors raise *predicted* height but not
reliably *attained* height; vosoritide raises velocity with final height still unknown. This
surrogate-unreliability is captured quantitatively in `surrogate_calibration` and is the single most
important cross-cutting caution in the warehouse.
