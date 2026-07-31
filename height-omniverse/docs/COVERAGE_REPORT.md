# Coverage report (NOT a verdict)

Reports what was searched and measured coverage. It does **not** answer "what is the best stack?" or
"what should someone take?" — the corpus has no recommendation layer.

## Databases / APIs searched this session
NCBI E-utilities (PubMed/PMC), Europe PMC REST, ClinicalTrials.gov, NCBI GEO / BioProject, HGNC
(complete gene set), Crossref/DOI, OMIM (via literature), FDA/EMA (labels/press where reachable).
Retrieval date **2026-07-31**. Full API/query logs are embedded in the four harvest-agent transcripts.

## Query families
- Interventions × {adult height, near-final height, growth velocity, predicted adult height, bone age,
  final-height RCT} across endocrine, CNP/NPR2, FGFR3, PTH1R, IGF, mechanical, nutrient/botanical.
- Human genetics: GIANT/height GWAS; monogenic stature allelic series (LoF/GoF/dup/del); estrogen/
  fusion boundary; overgrowth syndromes.
- Mechanical: PEMF/LIPUS/ESWT/loading/compression/distraction/vibration + exact waveform metadata.
- Omics: named anchor accessions + GEO discovery for postnatal growth-plate/epiphyseal datasets.

## Counts currently in the warehouse
| entity | count |
|---|---|
| gene universe (HGNC approved) | 45,030 |
| genes with recorded direction-of-effect | 21 |
| interventions (curated, verified) | 14 |
| studies (mechanical/nutrient/genetics) | 11 |
| evidence edges | 48 |
| surrogate→final-height calibration rows | 17 |
| growth-plate omics datasets catalogued | 11 (incl. zonal LMD + 2026 human) |
| access failures logged | 11 |
| evidence cards generated | 49 |

## Human final-height coverage (the scarce, decisive layer)
Interventions with **verified human attained-adult-height** evidence: GH (ISS PMID 15240584; GHD
10102057), oxandrolone+GH (Turner 21493672), GnRH agonist (early-onset CPP PMC5752503; normal-puberty
ISS 12621135), estrogen adjunct (Turner 21449786), mecasermin (SPIGFD 40626687, partial), levothyroxine
(11580992). **Surrogate-did-not-translate** documented for aromatase inhibitors in ISS boys (31024444)
and normal-puberty GnRHa. **Velocity-only / final-height UNKNOWN:** vosoritide, navepegritide,
infigratinib, TYRA-300, long-acting GH.

## Species/site/endpoint coverage
Heavy mouse bias in omics (only human postnatal GP dataset = GSE288028, raw withheld). Mechanical
final-length data exist essentially only in rat (LIPUS null; Hueter-Volkmann). No mechanical
intervention shows a **final adult bone-length** increase in humans.

## Saturation (honestly: NOT reached)
The spec's saturation criterion (3 discovery rounds with no new Tier-1/2 intervention, human
final-height study, growth-plate dataset, or dual-evidence gene) was **not run to closure** in this
session. One structured seed round + targeted recursive expansion was completed across four agents.
Measured coverage is a *starting* corpus, not a saturated one. Unresolved regions are logged in
`BLOCKED_DATA.md` and `LIMITATIONS.md`.

## Unresolved blind spots
LINCS/CMap perturbational layer (not ingested); genome-wide omics reprocessing from raw (only prior
GSE288529 audit reused); ChEMBL/PubChem/BindingDB compound-affinity layer (not ingested); patents/
theses/conference deep-mining (partial); restricted individual-level genetics (public summary only);
formal fine-mapping / MR / cell-type LDSC-SEG on GWAS (not run — GWAS captured at anchor level).
