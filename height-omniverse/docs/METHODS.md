# METHODS & hard rules

## Question & primary endpoint
Primary biological endpoint: **ΔL_stable_final = mature length (intervention) − mature length
(control)**, measured after washout, mineralization, remodeling, and skeletal maturity. Lower
endpoints (velocity, PAH, plate thickness, bone age, condyle, repair, expression) are stored but
**never** allowed to stand in for a higher one (`config/endpoint_ontology.yaml`).

## Research-only rule
Doses/PK/routes are captured as **scientific metadata only**. The corpus contains no human protocol,
dosing schedule, sourcing, injection, or DIY mechanical routine, and no "best stack" or recommendation
field anywhere.

## Evidence sourcing
- Conclusions rest on **primary** research / regulatory documents / registries / original datasets.
  Reviews were used only for discovery.
- Every load-bearing number traces to a **source_id** (PMID/PMCID/DOI/NCT/accession/URL). Values behind
  paywall/OCR are marked `NOT RETRIEVED`; memory-only statements are marked `UNVERIFIED(memory)` and
  are non-load-bearing.
- Retrieval via NCBI E-utilities, Europe PMC REST, ClinicalTrials.gov, GEO, HGNC, and publisher/
  registry pages. Retrieval date 2026-07-31.

## Directness tiers (`config/evidence_tiers.yaml`)
T1 human final height under exact intervention · T2 fusion-capable animal final length · T3 exact
intervention, conventional physis, no maturity · T4 related ligand/target in physis · T5 human
genetics · T6 condyle/fracture/organoid/non-fusing species · T7 cell/molecular markers · T8
simulation/inference/anecdote. Tiers are never naively pooled.

## Quality-control rules (enforced throughout)
- Cells are **never** biological replicates; one-library pseudobulk is descriptive only.
- mRNA ≠ protein ≠ pathway activity; dropout ≠ absence.
- No pooling of incompatible tissues (condyle ≠ tibia), of repair with growth, or of embryonic
  patterning with postnatal augmentation without a transportability model.
- Plate thickness / open plate / delayed ossification ≠ greater final length.
- Animal dose is never translated into a human recommendation.
- Nulls and inconvenient findings are preserved; malformed/weak bone is flagged when length rises.
- Species-fusion flag required before calling any rodent result "final-height evidence".
- Efficacy and safety analyzed independently, then connected only in a separate net-feasibility view.
- No single favored theory organizes the schema: fixed/declining reserve, population asymmetry,
  estrogen-accelerated senescence, local output-per-cell, geometry/matrix, niche plasticity,
  site-specific closure, and endocrine/mechanical spatial encoding are all representable.

## Surrogate calibration (the load-bearing analysis)
For interventions with both surrogate and final-height data we record whether the surrogate
**translated** to attained adult height (`surrogate_calibration` table). This directly quantifies the
`P(adult-height gain | velocity/PAH/BA change)` reliability that the endpoint hierarchy demands.

## Reproducibility
DuckDB warehouse rebuilt from `data/raw` + `data/curated` by `src/ingestion/build_db.py`; one-command
`run_all.sh`; pinned env (`environment/python.lock`); seed 20260731; raw checksums in
`docs/SOURCE_MANIFEST.csv`. Simulation outputs are always labeled simulation.

## Honest scope
See `docs/COVERAGE_REPORT.md` (what was searched + saturation), `docs/BLOCKED_DATA.md` (inaccessible
sources), and `docs/LIMITATIONS.md` (what was scaffolded but not run). No completeness is claimed.
