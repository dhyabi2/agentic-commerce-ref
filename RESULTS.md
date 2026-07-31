# HDAC4 repressive activity across growth plate zones and developmental age

**Hypothesis under test.** Growth plate closure reflects collapse of a self-reinforcing
repression state centred on nuclear HDAC4, rather than exhaustion of a progenitor budget.
The prediction is that HDAC4 repressive activity in resting and proliferative chondrocytes
declines with developmental age, and that the zone boundary where repression is lost
migrates toward the resting zone.

**Critical limitation, stated up front and repeated throughout.** HDAC4 is regulated by
nuclear/cytoplasmic shuttling — a post-translational process. RNA-seq cannot measure
subcellular localisation. Everything below infers repressive activity from the expression
of HDAC4 target genes. A negative result does not refute the hypothesis; a positive result
is consistent with, but does not demonstrate, changed localisation. The definitive
experiment is immunostaining for HDAC4 subcellular distribution across a developmental time
course in a fusing species.

---

## The finding

**The temporal prediction cannot be tested in the specified human data, and the reason is
not low power — it is that the regressor does not exist.**

GEO carries **no per-donor age, sex, or Tanner stage** for GSE288028. The sample records
contain exactly five characteristic tags — `tissue`, `cell type`, `genotype`, `treatment`,
`batch` — verified in both the SOFT sample records and the MINiML family XML. The only age
information published is one sentence in the series-level overall design: *"four human
adolescent growth plate chondrocytes were isolated from patients with age 11-14 years."*
That is a group range with no mapping to individual donors, so it cannot be joined to any
sample. Analysis B was therefore **not fitted**, per the pre-registered adequacy gate. No
surrogate — donor index, accession order, batch order — was substituted.

The second named dataset, GSE288529, is **one animal at one age** (female C57BL/6J, 4
weeks). The two mouse libraries inside GSE288028 are two animals of unrecorded age from a
single batch. Neither contributes an age contrast.

A systematic GEO search (18 logged queries, 462 distinct series, 146 growth-plate scRNA-seq
series, `results/tables/geo_search_log.tsv`) found **no growth plate scRNA-seq dataset with
a developmental age series** — not for human, not for mouse. The nearest misses are
enumerated in section 1.

Even if per-donor ages were recovered from the authors, they could not be interpreted:
donor is perfectly confounded with processing batch **and** with CellRanger version (6.0.1
for D3/D4, 7.1.0 for D1/D2), so an age effect and a pipeline-version effect are not
separable in this design.

What *can* be said is set out below. In brief: the activity score behaves as a working
score should (section 2); the only interventional contrast in the data is null (section 5);
and the only dataset in GEO with genuine zone × age structure — rat, bulk microarray, from
2009 — shows target expression rising with age in the direction the hypothesis predicts,
but that result rests on its oldest age group alone (section 3.2).

---

## 1. What datasets exist, and what age span they cover

| Dataset | Organism | Design | Age information | Usable age contrast |
|---|---|---|---|---|
| **GSE288028** (primary) | human + mouse | 4 donors; uncultured / 24 h vehicle / 24 h GH; 12 human + 2 mouse libraries | **none per donor**; series-level "11–14 years" only | **No** |
| **GSE288529** | mouse | 1 animal, epiphysis digest | 4 weeks, single value | **No** |
| GSE16981 (supplementary) | rat | microdissected zones, n=5 per group | **1, 3, 6, 9, 12 weeks, per sample** | **Yes** — see 3.2 |

Post-QC cell counts: human 67,643 cells (12 libraries); mouse GSE288028 11,352 (2);
mouse GSE288529 11,816 (1). Full QC in `results/tables/qc_per_sample.tsv`.

### 1.1 Search for a growth plate scRNA-seq age series — result: none exists

18 queries logged with hit counts in `results/tables/geo_search_log.tsv`; all 462 distinct
series in `geo_search_hits.tsv`; the 146 growth-plate scRNA-seq series in
`geo_search_candidates.tsv`. The closest candidates, each checked individually:

- **GSE76157** — "growth plate development ... single-cell RNA-seq", 217 cells. All from
  **postnatal day 7**. One age.
- **GSE244880/GSE244881** — resting zone scRNA-seq, murine femur, **P36**. One age.
- **GSE267139** — "Human cartilage development at single cell resolution", 25 samples.
  hPSC-derived chondrocytes plus **one fetal femur at E59**. Not a postnatal age series.
- **GSE182540** — growth plate/periosteum stem cells in a *Zmpste24* premature-ageing
  **genetic model**, not a chronological age series.
- **GSE18338** — human growth plate at **Tanner B2 and B3 in one girl**. A genuine
  within-subject pubertal progression, but bulk microarray, n=1 patient.
- **GSE16981** — rat, zones microdissected, **1/3/6/9/12 wk**, n=5. Bulk microarray. The
  only public dataset with real zone × age structure; analysed in 3.2.

**Is the available age span adequate to address the question?** For the specified human
data: no — there is no span, because there are no ages. For the field as a whole: no public
single-cell dataset yet spans human growth plate development in a way that could test this
prediction.

### 1.2 The human libraries barely contain hypertrophic chondrocytes

| Uncultured library | Hypertrophic cells recovered |
|---|---|
| D1_P30453 | 11 |
| D2_P31011 | 3 |
| D3_P25452 | 9 |
| D4_P22202 | 1 |

20 of 60 human zone × library combinations fall below the pre-registered 100-cell floor
(0 of 15 for the mouse datasets). D4_P22202 yields only 375 cells in total, so **every** one
of its zones is flagged. Applying the pre-registered filter strictly leaves a single human
donor with three usable zones, which is why section 2 reports the filtered and unfiltered
passes side by side. Hypertrophic chondrocytes are large and fragile and are known to be
lost in droplet dissociation; this is a property of the assay, not of these donors.

---

## 2. Zonal sanity check — the score works

Zones were assigned from marker genes only (`SFRP5/APOE/CLU/PTCH1/CYTL1/GAS1/PTHLH/RAMP3`;
`CCND1/MKI67/E2F1`; `IHH/ALPL`; `COL10A1/MMP13`), independent of the scoring genes. Three
genes appear in both lists — **ALPL, COL10A1, MMP13** — and all three were excluded from the
primary score before fitting, leaving `IBSP, SPP1, VEGFA, PANX3, SP7`. (The brief flags
"IHH and COL10A1" as the overlap; IHH is in fact not in the target set, and ALPL and MMP13
are. All actual overlaps were removed.)

Mean pseudobulk target expression by zone (primary set):

| Dataset | resting | proliferative | prehypertrophic | hypertrophic |
|---|---|---|---|---|
| human GSE288028 | −0.039 | −0.125 | **0.566** | 0.510 |
| mouse GSE288028 | 0.011 | −0.011 | **0.366** | 0.196 |
| mouse GSE288529 | −0.023 | −0.224 | 0.406 | **0.563** |

Rank correlation with zone order, aggregated to the donor (the biological replicate):

| Dataset | method | mean ρ across donors | 95% CI | n donors | donors positive |
|---|---|---|---|---|---|
| human | scanpy `score_genes` | **0.690** | 0.317 – 1.063 | 4 | 4/4 |
| human | decoupler ULM | **0.727** | 0.307 – 1.146 | 4 | 4/4 |
| mouse GSE288028 | both | 0.600 | — | 2 | 2/2 |
| mouse GSE288529 | scanpy / ULM | 0.800 / 0.600 | — | 1 | 1/1 |

**The two methods agree in direction everywhere.** No disagreement to report.

**Verdict: the pre-registered gate passes, with one honest qualification.** Both operational
criteria are met — ρ > 0, and hypertrophic > resting in all three datasets. But the gradient
is **not monotonic**: in human and in mouse GSE288028 it peaks at the *prehypertrophic* zone
and falls slightly at the terminal hypertrophic zone (see `fig1_zonal_gradient.png`). The
inverted-U is expected here — the primary set is dominated by prehypertrophic/osteogenic
genes (IBSP, SPP1, SP7) after the three terminal markers were removed for non-circularity,
and the human hypertrophic pseudobulks rest on 1–11 cells each. The claim supported is
"target expression rises across the differentiation axis and is higher in hypertrophic than
in resting cells", not "target expression increases monotonically".

Because the gate passes, downstream analyses are interpretable — subject to everything in
section 7.

### 2.1 Negative controls (all pre-registered)

| Control | Result | Pass? |
|---|---|---|
| Housekeeping set, same procedure | mean ρ = −0.167 (scanpy) / −0.033 (ULM); permutation p = 0.32 / 0.87 — inside the null band | **Yes** |
| Zone-label shuffle, 1000 permutations | observed ρ 0.617 (scanpy) / 0.667 (ULM); null mean −0.049, null 95% band [−0.35, 0.27]; **empirical p = 0.001** (the floor at 1000 permutations) | **Yes — gradient vanishes** |
| Per-zone cell counts reported, <100 flagged | 20/60 human, 0/15 mouse combinations flagged; all reported | **Yes** |
| Leave-one-donor-out | primary ρ 0.578 – 0.667 (scanpy), 0.600 – 0.756 (ULM) across all four leave-outs; housekeeping stays −0.20 to +0.06 | **Yes — no single donor drives it** |
| Score vs library depth / mitochondrial fraction | primary score: all \|ρ\| ≤ 0.34, p ≥ 0.069 | **Yes** |

One control result deserves emphasis rather than a tick. The **housekeeping** score *is*
strongly correlated with technical covariates — total counts ρ = 0.46 (p = 0.012), n_genes
ρ = 0.40 (p = 0.030), % mitochondrial ρ = 0.54 (p = 0.0023) — while the primary score is
not. The housekeeping set is behaving exactly as a technical-artefact detector should, and
it says the primary score is comparatively clean. But the housekeeping profile is not
perfectly flat either (1.01 → 1.15 → 1.26 → 1.13 across zones, `fig5_controls.png`); its
range is roughly a third of the primary score's and its rank correlation is
indistinguishable from the permutation null, but some small part of the zonal profile is a
general shift in library complexity across zones, not HDAC4 target biology.

---

## 3. Age effect

### 3.1 Human — UNDEFINED, not fitted

Pre-registered gate (a) — per-donor age recorded — **fails**. No model was fitted.
`results/tables/B_age_adequacy_gate.json` records the verdict machine-readably.

What can be quantified is the donor-to-donor spread any age effect would have had to exceed.
Resting-zone pseudobulk, uncultured libraries, primary set, scanpy scoring:

- n = 3 donors passing the 100-cell floor (D4 excluded at 44 resting cells)
- mean −0.175, **between-donor SD 0.111**, range 0.206 (−0.255 to −0.048)
- 95% CI half-width on the donor mean at this n: **±0.275**
- minimum detectable difference at 80% power, α = 0.05, n = 3: **1.27** — roughly **6×**
  the entire observed spread between donors

For contrast, the **technical** interval from bootstrapping cells within a single library
(1000 resamples, labelled `TECHNICAL` throughout `B_resting_zone_TECHNICAL_bootstrap.tsv`):

| Donor | resting cells | mean | TECHNICAL 95% CI | width |
|---|---|---|---|---|
| D1_P30453 | 1292 | −0.221 | −0.232 to −0.209 | 0.023 |
| D2_P31011 | 694 | −0.255 | −0.274 to −0.234 | 0.040 |
| D3_P25452 | 4072 | −0.048 | −0.057 to −0.039 | 0.018 |
| D4_P22202 | 44 | −0.235 | −0.298 to −0.159 | 0.139 |

The technical intervals are **~5–10× narrower than the between-donor SD**. Treating cells as
replicates here would have produced intervals that exclude zero for almost any comparison
while the actual donor-level evidence is nil. This is the concrete reason the pre-registration
forbids it.

### 3.2 Rat GSE16981 (supplementary) — the only real age series

Not part of the primary analysis. Rat, bulk Affymetrix (MAS5), zones separated by manual
microdissection, n = 5 biological replicates per group. Animals older than 1 week were
**castrated at day 18**; the 3–12 wk series is internally consistent (all castrated, same
combined proliferative + early-hypertrophic dissection), so the age slope within it is not
confounded by gonadal status — but these are gonadal-steroid-deficient animals throughout.
The 1 wk arm cannot be pooled in: different dissection (separate zones) and intact gonads.

**Zonal sanity check at 1 wk** (all three zones present): ρ = 0.661 (p = 0.0072) by z-mean,
ρ = 0.832 (p = 0.00012) by ULM; housekeeping ρ = −0.227 (p = 0.42). The score behaves the
same way in rat bulk data as in the single-cell data.

**Age effect, proliferative + early hypertrophic zone, castrated animals, 3 → 12 weeks:**

| Score | slope / week | 95% CI | p | R² | n arrays |
|---|---|---|---|---|---|
| primary, z-mean | **+0.0454** | +0.0225 – +0.0682 | **0.00059** | 0.49 | 20 |
| primary, decoupler ULM | +0.0374 | +0.0114 – +0.0634 | 0.0073 | 0.34 | 20 |
| secondary (8-gene), z-mean | +0.0298 | +0.0074 – +0.0521 | 0.012 | 0.30 | 20 |
| **housekeeping control** | −0.0132 | −0.0463 – +0.0199 | 0.41 | 0.04 | 20 |

Positive slope = target expression rises with age = **inferred HDAC4 repressive activity
falls with age**. That is the direction the hypothesis predicts, in the proliferative
compartment, with a flat housekeeping control and agreement between both scoring methods.

**And here is the inconvenient part.** Leave-one-age-out shows the slope is carried by the
oldest group:

| Age left out | slope / week | p |
|---|---|---|
| none (full) | +0.045 | 0.00059 |
| 3 wk | +0.045 | 0.032 |
| 6 wk | +0.049 | 0.00098 |
| 9 wk | +0.053 | 0.000055 |
| **12 wk** | **+0.020** | **0.245** |

Group means are −0.203 (3 wk), +0.011 (6 wk), −0.083 (9 wk), +0.282 (12 wk) — not a smooth
progression. Drop the 12 wk animals and there is **no significant trend across 3–9 weeks**.
The honest summary is: a positive age slope over 3–12 weeks in castrated rat, driven by the
oldest age group, in a compartment that is a mixture of proliferative and early hypertrophic
cells, measured as bulk transcript abundance. It is consistent with the prediction. It is
one age group away from being nothing, it is not human, it is not single-cell, and rats do
not undergo the human fusion event.

---

## 4. Boundary migration

Per-library differentiation pseudotime (diffusion pseudotime rooted at the highest
resting-marker cell), with the crossing point defined as the position where the target score
first exceeds the midpoint between that library's resting and hypertrophic means.

**Cannot be tested against age** — same missing regressor as section 3.1. But the crossing
positions were computed and reported per library anyway, and they carry their own negative
result:

| Donor | arm | crossing quantile |
|---|---|---|
| D1_P30453 | uncultured | 0.67 |
| D2_P31011 | uncultured | 0.07 |
| D3_P25452 | uncultured | 0.41 |
| D3_P25452 | vehicle (file 1) | 0.03 |
| D3_P25452 | vehicle (file 2) | 0.29 |
| D3_P25452 | GH (file 1) | 0.73 |
| D3_P25452 | GH (file 2) | no crossing |

**The metric is not stable enough to detect anything.** The same donor, same arm, split
across two sequencing libraries gives 0.03 versus 0.29 (vehicle) and 0.73 versus no crossing
at all (GH). Two of thirteen libraries have no crossing because their hypertrophic mean sits
below their resting mean. The crossing quantile is a *rank* within each library, so it is
dominated by that library's zone composition — which, as section 1.2 shows, varies by an
order of magnitude between libraries of the same donor. Even with per-donor ages in hand,
this metric as specified could not have supported a conclusion in this dataset.

---

## 5. GH contrast — the only interventional comparison, and it is null

GH versus vehicle, paired within donor, resting/progenitor cells, n = 3 donors (D4 has no
culture arm):

| Gene set | Method | Restricted to GP2 | GH − vehicle (target) | 95% CI | p |
|---|---|---|---|---|---|
| primary | scanpy | yes | +0.076 | −0.549 – +0.701 | 0.65 |
| primary | ULM | yes | +0.169 | −2.928 – +3.266 | 0.84 |
| primary | scanpy | no | +0.034 | −0.426 – +0.495 | 0.78 |
| housekeeping | scanpy | yes | −0.051 | −0.165 – +0.062 | 0.19 |

Per-donor differences run in opposite directions (D1 +0.35, D2 +0.03, D3 −0.15). **No
detectable shift in the HDAC4 target score under GH.** The confidence interval is wide
enough to contain effects several times larger than the entire between-zone range, so this
is uninformative rather than evidence of absence. It is **not** a trend.

**The pre-specified GP1 confound is worse than stated in the brief.** The brief warns that
24 h culture loses the GP1 quiescent cluster, so the contrast should be restricted to GP2.
Subclustering the resting-zone cells shows that culture does not remove one subcluster — it
separates *all* of them (`D_rz_subcluster_culture_depletion.tsv`):

| RZ subcluster | uncultured cells | cultured cells |
|---|---|---|
| 3 | 1061 | 3 |
| 2 | 3003 | 26 |
| 0 | 1801 | 26 |
| 4 | 232 | 312 |
| 1 | 3 | 620 |
| 6 | 1 | 1887 |
| 7 | 1 | 5650 |
| 5 | 0 | 650 |
| 8 | 0 | 1244 |

Only subcluster 4 contains meaningful numbers of both. Twenty-four hours of organ culture
reshapes the resting-zone transcriptome so completely that cultured and uncultured cells
share essentially no subcluster. "Restrict to GP2" is therefore not an available fix — every
uncultured subcluster is lost, not just the quiescent one. The GH-versus-vehicle contrast
itself remains valid because **both** arms are cultured and share that baseline shift, which
is why the restriction changes the answer barely at all (+0.076 versus +0.034). But no
statement about the uncultured resting zone can be carried over to the cultured arms.

---

## 6. Pathway-component expression by zone

Reported directly; never folded into the activity score. Human GSE288028, mean
log-normalised expression:

| Gene | resting | proliferative | prehypertrophic | hypertrophic |
|---|---|---|---|---|
| **HDAC4** | **0.417** | 0.309 | 0.263 | **0.153** |
| HDAC5 | 0.195 | 0.221 | 0.174 | 0.147 |
| HDAC7 | 0.100 | 0.132 | 0.100 | 0.072 |
| SIK1 | 0.001 | 0.001 | 0.000 | 0.000 |
| SIK2 | 0.438 | 0.517 | 0.169 | 0.262 |
| SIK3 | 1.006 | 0.960 | 0.717 | 0.846 |
| CAMK4 | 0.049 | 0.117 | 0.033 | 0.074 |
| CAMK2D | 0.638 | 0.643 | 0.507 | 0.484 |
| PPP2CA | 0.377 | 0.488 | 0.289 | 0.222 |
| PPP2R1A | 0.362 | 0.411 | 0.313 | 0.242 |
| PTH1R | 0.916 | 0.547 | 1.125 | 1.136 |
| PTHLH | 0.047 | 0.031 | 0.006 | 0.003 |
| **MEF2C** | 0.495 | 0.675 | **0.936** | 0.708 |
| **RUNX2** | 0.895 | 0.447 | **1.343** | 0.944 |
| SMAD2 | 0.499 | 0.526 | 0.357 | 0.451 |
| SMAD3 | 0.273 | 0.280 | 0.315 | 0.365 |
| HIF1A | 1.413 | 1.246 | 0.879 | 1.010 |
| EPAS1 | 0.334 | 0.667 | 0.716 | 0.380 |

Descriptively coherent with the model: **HDAC4 transcript is highest in the resting zone and
declines monotonically to the hypertrophic zone** (0.417 → 0.153, a 2.7-fold drop), PTHLH
falls 16-fold along the same axis, and MEF2C and RUNX2 both peak in the prehypertrophic zone
where the target genes come on. The same pattern holds in both mouse datasets
(`fig4_pathway_heatmap.png`).

This is suggestive and it is also exactly the observation the central limitation applies to.
Falling *HDAC4 transcript* toward the hypertrophic zone is not falling *nuclear HDAC4*. A
cell could hold HDAC4 mRNA constant and de-repress its targets entirely by exporting the
protein — that is the mechanism the hypothesis actually proposes, and none of these numbers
can see it.

---

## 7. Limitations

1. **Transcript ≠ protein ≠ localisation.** The whole analysis is an inference from target
   transcript abundance. HDAC4 acts by nuclear/cytoplasmic shuttling. Nothing here measures
   or can measure that. This is not a caveat on the result; it is a bound on what the assay
   is capable of.
2. **No age variable.** The primary comparison was undefined, not underpowered. Nothing in
   this report should be read as a test of the temporal prediction in human tissue.
3. **Donor n = 4, and 3 after the cell-count floor.** The sign test cannot reach p < 0.05
   with four donors even at 4/4 concordance (minimum two-sided p = 0.125). Absence of
   significance in this dataset carries almost no evidential weight.
4. **Donor is confounded with batch and with CellRanger version** (6.0.1 vs 7.1.0). Any
   between-donor contrast, age or otherwise, is confounded with pipeline version.
5. **The human hypertrophic zone is effectively absent** — 1 to 11 cells per uncultured
   library. The upper end of the differentiation axis is carried by a handful of cells.
6. **Mouse ≠ human, and mice do not fuse.** The mouse datasets can validate the score; they
   cannot inform a claim about growth plate closure. No mouse result is extrapolated to
   human fusion anywhere in this report.
7. **Rat ≠ human either**, and the rat animals were castrated. The supplementary age result
   comes from gonadal-steroid-deficient rats, in bulk, in a mixed proliferative/early-
   hypertrophic compartment, and it depends on its oldest age group.
8. **Pubertal stage range.** Even taking the series-level statement at face value, all four
   human donors sit inside a single narrow pubertal window. Had per-donor ages existed, that
   window would probably still have been too narrow to test a prediction about closure.
9. **Culture artefacts are severe.** 24 h organ culture separates essentially every
   resting-zone subcluster from its uncultured counterpart. Cultured and uncultured arms are
   not interchangeable for any purpose.
10. **Raw reads for the human samples were never deposited** (withheld for participant
    privacy). This analysis inherits the submitters' alignment and cell-calling decisions
    and cannot audit them.
11. **The primary score is 5 genes.** Removing ALPL, COL10A1 and MMP13 to avoid circularity
    with zone assignment was the right call, but it leaves a small set weighted toward
    prehypertrophic/osteogenic biology, which is why the zonal profile peaks before the
    terminal zone.
12. **No cross-library batch integration was applied.** Highly variable genes were selected
    with `batch_key`, but no Harmony/BBKNN/scVI correction was run, so the human embedding
    separates strongly by library (`fig6_umap_zones.png`, 40 leiden clusters across 12
    libraries) and clusters partly track library of origin. Zone assignment is
    marker-driven per cluster and pseudobulk inference is done per library, so this does not
    invalidate the zonal result — the leave-one-donor-out and permutation controls both hold
    — but cluster identity should not be read as a shared cell state across donors, and the
    UMAP should not be read as an integrated atlas.

---

## 8. What would need to be measured to actually answer the question

1. **The definitive experiment: HDAC4 immunostaining for subcellular distribution**, scored
   as nuclear/cytoplasmic ratio per cell, per zone, across a developmental time course in a
   **fusing species** — human epiphysiodesis specimens across Tanner stages, or a species
   that closes its plates. No amount of RNA-seq substitutes for this. Everything in this
   report is a proxy for the one measurement that would settle it.
2. **Human growth plate scRNA-seq with per-donor age and Tanner stage actually recorded**,
   spanning pre-pubertal through post-fusion. The tissue is obtainable — epiphysiodesis is
   performed across a wide age range. The binding constraint is that existing studies did
   not deposit ages, not that the samples are impossible to get. Depositing age as a
   sample characteristic costs nothing and would have made this analysis possible.
3. **Enough donors to see past donor variance.** Between-donor SD in the resting zone is
   0.111 on the primary score. Detecting a difference of that size (1.0 SD) at 80% power,
   α = 0.05, needs **≈17 donors per group** (two-sample t, `statsmodels`: 16.7) — not 4.
4. **Spatial transcriptomics or microdissected zones**, so zone assignment does not depend on
   marker-based clustering of dissociated cells, and so hypertrophic chondrocytes are not
   lost to dissociation. GSE16981's microdissection approach applied to human tissue with a
   real age series would answer the zonal half of the question directly.
5. **Paired protein and transcript from the same specimens**, to establish whether target
   transcript abundance tracks HDAC4 nuclear localisation at all. Without that calibration,
   the sign of a target-expression result cannot be converted into a statement about
   repressive activity — which is the assumption the entire scoring approach rests on.
6. **A conditional nuclear-anchored HDAC4 allele** across the same time course, to test
   whether forcing HDAC4 nuclear retention delays closure. That distinguishes the
   repression-collapse model from the progenitor-exhaustion model causally, which no
   observational transcriptomic comparison can do.

---

## Reproducing

```bash
./run_all.sh          # fetches ~760 MB from GEO, runs everything, writes results/
```

Environment pinned in `requirements.txt`; seeds fixed in `src/hdac4_config.py`;
download dates, accessions and sha256 checksums in `docs/DATA_PROVENANCE.md`;
pre-registration in `docs/PREREGISTRATION.md`; departures from it in `docs/DEVIATIONS.md`.

| Figure | Content |
|---|---|
| `fig1_zonal_gradient.png` | zonal score gradient per dataset, both methods, per-library lines |
| `fig2_score_vs_age.png` | human (no age axis exists) and rat (real age series) |
| `fig3_boundary_crossing.png` | boundary crossing position per library |
| `fig4_pathway_heatmap.png` | pathway-component expression by zone |
| `fig5_controls.png` | housekeeping control and permutation null |
| `fig6_umap_zones.png` | zone assignment |
