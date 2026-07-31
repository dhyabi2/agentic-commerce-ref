# Data dictionary

DuckDB warehouse `database/height_omniverse.duckdb`. Rebuild: `src/ingestion/build_db.py`.

## `genes` (genome-wide universe; 45,030 approved HGNC genes)
`ensembl_id, hgnc_id, symbol, name, locus_group, locus_type, biotype, chromosome, aliases,
prev_symbols, gene_group, entrez_id, status` (from HGNC complete set) +
`height_gwas_evidence, monogenic_stature_evidence, growthplate_zone_evidence, direction_of_effect,
source_provenance` (curated; sparse — populated only where evidence exists). Mirror: `data/processed/ALL_GENES.parquet`.

## `interventions`
`intervention_id, name, synonyms, modality, target, action, approval_status, route, pk_summary,
evidence_context, technology_readiness, source_provenance`. Modality ∈ {hormone, small molecule,
peptide, antibody, fusion protein, oligonucleotide, gene therapy, editing, cell therapy,
tissue-engineered, nutrient, botanical, mechanical waveform, device, exercise, surgery,
disease correction, environmental}.

## `studies`
`study_id, citation, year, design, species, sex, age_stage, disease, tissue_site, intervention,
comparator, n_true_bioreplicates, exposure, duration, washout, endpoint, maturity_followup,
effect_size, uncertainty, direction, geometry_mineralization, directness_tier, risk_of_bias,
source_location, notes`. `directness_tier` ∈ T1..T8 (`config/evidence_tiers.yaml`).

## `outcomes`
`outcome_id, study_id, outcome_type, raw_value, standardized_value, unit, endpoint_rank,
is_final_length, notes`. `endpoint_rank` from `config/endpoint_ontology.yaml` (1 = final length).

## `surrogate_calibration` (load-bearing)
`intervention, population, design, n, surrogate_effect, adult_height_effect, surrogate_translated,
direction, source_id, notes`. `surrogate_translated` ∈ {YES, NO, UNKNOWN, PARTIAL} — did a
velocity/PAH/bone-age gain become ATTAINED adult height?

## `evidence_edges` (knowledge graph)
`edge_id, subject, predicate, object, direction, context, effect_size, source_id, evidence_tier,
directness, confidence, conflicting_edges`. Exported to `knowledge_graph/*.graphml|parquet`.

## `evidence_snippets`
`snippet_id, source_id, excerpt, paraphrase, location, extractor_confidence, human_reviewed, links_to`.

## `access_failures`
`source, reason, attempted_method, metadata_found, request_route, importance`. Every blocked source is
recorded here and in `docs/BLOCKED_DATA.md` — no silent source failure.

## Provenance rule
Every load-bearing value carries a `source_id` (PMID/PMCID/DOI/NCT/accession/URL). `NOT RETRIEVED` =
value behind paywall/OCR; `UNVERIFIED(memory)` = non-load-bearing, not used for conclusions.
