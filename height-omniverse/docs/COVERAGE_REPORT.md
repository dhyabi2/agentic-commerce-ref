# Coverage & saturation audit (measured — NOT a verdict)

Reports what was searched and measured coverage with **per-category saturation labels**. It does not
answer "what is the best stack?" — the corpus has no recommendation layer.

Label key: **saturated** (3 rounds, no new qualifying item) · **not saturated** (real items remain) ·
**externally blocked** (public data insufficient/withheld) · **structurally unknowable** (cannot be
resolved from public data).

## Databases / APIs actually queried (2026-07-31)
NCBI E-utilities (esearch/esummary/efetch/elink), Europe PMC REST, ChEMBL REST, NCBI GEO/BioProject,
HGNC complete set, Crossref/DOI, OMIM (via literature). Query/round logs:
`results/tables/SEARCH_QUERIES.csv`, `SCREENING_LOG.csv`, `CITATION_FRONTIER.csv`, `PRISMA_COUNTS.csv`.

## Warehouse contents now
| entity | count | note |
|---|---:|---|
| gene universe (HGNC approved) | 45,030 | genome-wide |
| genes with direction-of-effect | 21 | verified PMIDs |
| genome-wide zonal atlas rows | 98,788 | 24,697 genes × 4 states (GSE288529) |
| interventions | 14 | verified |
| literature studies (screened) | 274 incl | 4 citation rounds; PRISMA ~1200 identified |
| surrogate→final-height calibration | 17 | crown jewel |
| compound–target activities | 275 fetched / 19 named | ChEMBL, 10 targets w/ data |
| evidence edges | 48 | knowledge graph |
| interaction pairs | 91 (14 non-trivial) | pathway convergence |
| dual-evidence high-priority genes | 11 | ≥3 modes |
| QA tests | 40 passing | +saturation acceptance |

## Per-category saturation

| category | status | basis |
|---|---|---|
| Schema/QA integrity | **saturated** | 34 tests; migration hardened; all edges/studies source-linked |
| Human final-height interventions (endocrine) | **partially saturated** | 4 citation rounds run; **274 included studies**; 5/15 categories saturated; clinical/final-height categories converged; Tier-7/8 mechanism categories (NPR2/CNP, growth-plate biology, ACAN, limb-lengthening) capped at 4 rounds with open residual logged (docs/SATURATION_ANALYSIS.md) |
| Surrogate→final-height calibration | **not saturated** (rich) | 17 curated; velocity-only class (vosoritide/CNP/FGFR3i) growing |
| Human height genetics (monogenic/allelic) | **not saturated** | 21 verified genes + GWAS anchors; full GWAS-sumstat fine-mapping not run |
| GWAS common-variant fine-mapping / MR / LDSC-SEG | **externally blocked (compute/data)** | needs LD panels + full sumstats; enrichment done on curated set (perm p=0.023) |
| Growth-plate omics (mouse) | **not saturated** | GSE288529 reprocessed genome-wide; 10+ more catalogued, not reprocessed from raw |
| Growth-plate omics (human postnatal) | **externally blocked** | only GSE288028; human raw withheld (PII) |
| Compound–target (ChEMBL) | **not saturated** | 10 targets w/ data ingested; PAPPA2/GHR/NPR2/IGF1 have no small-molecule bioactivity (structural) |
| LINCS/CMap perturbational | **externally blocked (size)** | ~50 GB L1000; reproducible fetch script stubbed; not ingested |
| Regulatory/patent/grey full-text | **not saturated / partially blocked** | abstracts captured; FDA/EMA PDFs not text-extractable this pass |
| Spatial reaction-diffusion model | **saturated (model class)** | runs; mass-conservation + gradient-recovery tests pass |
| Agent-based lineage model | **built** | neutral drift -> monoclonality (Newton 2019); reserve-depletion -> fusion |
| Finite-element mechanical model | **externally blocked** | FEniCSx not installed; spatial RD + Hueter-Volkmann evidence serve as surrogate (backlog J1) |
| Mechanical final-length (human) | **structurally unknowable (currently)** | no human mechanical final-length datum exists (LIPUS null; PEMF=thickness) |

## Cross-cutting measured signal (data, not advice)
- **Surrogate unreliability:** aromatase inhibitors raise predicted height but not attained adult
  height (PMID 31024444); GnRHa in normal-puberty ISS ≈ neutral; vosoritide/CNP/FGFR3i are
  velocity-only with final height unknown. GH, oxandrolone+GH, GnRHa-in-early-CPP, mecasermin, and
  levothyroxine translate.
- **Height genes are growth-plate-expressed** more than random (perm p=0.023), split resting
  (PAPPA2/IGF1/ESR1/NPR2) and prehypertrophic (FGFR3/IHH/ACAN); no single-zone enrichment reaches
  significance in the verified set — the Yengo-2022 resting-zone-enrichment claim is **not** overstated.
- **Systemic delivery flattens spatial information** and penetrates the resting zone poorly (spatial
  RD median penetration 0.29) — consistent across PK/PD, growth-plate, and spatial models.

## Not exhaustive
Saturation is **not** reached across most discovery categories (see labels). This is a substantially
expanded, queryable corpus, not a complete one. Continuation in `docs/EXACT_CONTINUATION_PROMPT.md`.
