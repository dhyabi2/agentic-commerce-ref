# Completion state

- **Branch / PR:** claude/height-omniverse-data-harvest / PR #2
- **Completed phases:** A (schema+QA), D (genome-wide zonal atlas), E (GWAS×zone enrichment),
  F (ChEMBL compound–target), G (dual-evidence discovery), H (spatial reaction-diffusion),
  K (interaction matrix), M (query views), O (coverage audit).
- **Partial:** B (1 seed round + 1 citation round; not saturated), C (abstracts only), N (cards regenerated, not yet enriched), D2 (only GSE288529 reprocessed), I (ODE+spatial cover it; no dedicated ABM).
- **Externally blocked:** F2 LINCS (~50 GB), J FEniCS (uninstalled), human omics raw (PII), GWAS fine-mapping (LD panels).
- **Next phase:** see docs/EXACT_CONTINUATION_PROMPT.md (finish B saturation first).
- **Tests:** 34 passing.
- **Datasets downloaded:** HGNC complete set (checksummed); ChEMBL activities (275 rows); GSE288529 processed h5ad reused (1.4 GB, external to repo).
- **Storage:** DB rebuildable (gitignored); parquet atlas ~ small; no large raw committed.
- **Unsuccessful/withheld:** GSE288028 human raw (PII); FDA/EMA PDFs (not text-extractable); PAPPA2/GHR/NPR2/IGF1 no ChEMBL small-molecule bioactivity.
