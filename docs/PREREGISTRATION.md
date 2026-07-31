# Pre-registration: HDAC4 repressive activity across growth plate zones and developmental age

**Status:** written before any score was computed and before any model was fit.
**Written:** 2026-07-31 (session date).

## 0. What was known when this document was written

This document was written **after** inspecting dataset *metadata* (GEO series records,
sample characteristic tags, matrix dimensions, CellRanger versions) and **before**
computing any activity score, fitting any model, or looking at any expression value.
That ordering is deliberate and is stated here so the reader can calibrate: the
adequacy gate in section 5 was written with knowledge of what age metadata does and
does not exist, because that knowledge is a property of the study design, not a result.

## 1. Hypothesis and its limits

**Hypothesis under test.** Growth plate closure reflects collapse of a self-reinforcing
repression state centred on nuclear HDAC4, rather than exhaustion of a progenitor budget.

**Prediction.** HDAC4 repressive activity in resting and proliferative chondrocytes
declines with developmental age, and the zone boundary at which repression is lost
migrates toward the resting zone.

**Critical limitation, restated in every output.** HDAC4 is regulated by
nuclear/cytoplasmic shuttling — a post-translational process. RNA-seq cannot measure
subcellular localisation. This analysis can only *infer* repressive activity from the
expression of HDAC4 target genes. A negative result does not refute the hypothesis, and
a positive result is consistent with, but does not demonstrate, changed localisation.
The definitive experiment is immunostaining for HDAC4 subcellular distribution across a
developmental time course in a fusing species.

## 2. Replicate definition (binding)

The **donor** (human) or **animal** (mouse) is the biological replicate. **Cells are
never treated as biological replicates.** All primary inference uses pseudobulk
aggregation to one value per sample per zone.

Any interval computed by bootstrapping or testing across cells within one library is a
**technical interval**. Wherever such an interval appears in a table, figure, or
sentence it is labelled `TECHNICAL` explicitly. It quantifies sampling noise within one
library and carries no information about donor-to-donor variability.

## 3. Zone-assignment markers (fixed)

Zone assignment must be independent of the scoring genes.

| Zone | Markers |
|---|---|
| resting / progenitor | SFRP5, APOE, CLU, PTCH1, CYTL1, GAS1, PTHLH, RAMP3 |
| proliferative | CCND1, MKI67, E2F1 |
| prehypertrophic | IHH, ALPL |
| hypertrophic | COL10A1, MMP13 |

### 3.1 Overlap between zone markers and scoring genes — resolution

The task brief flags IHH and COL10A1 as appearing in both lists. On direct comparison
of the two lists as specified, IHH is **not** in the scoring set; the genes that
genuinely appear in both are:

- **ALPL** (prehypertrophic marker ∩ repression-target set)
- **COL10A1** (hypertrophic marker ∩ repression-target set)
- **MMP13** (hypertrophic marker ∩ repression-target set)

All three overlapping genes are **excluded from the primary score** and retained for
zone assignment, per the brief's stated preference. This is recorded here, before
fitting, as a fixed choice.

Consequence, stated in advance: excluding the three strongest hypertrophic genes removes
most of the circularity from the zonal sanity check but also removes most of its
expected effect size. A weaker zonal gradient in the primary score than in the secondary
score is therefore an *expected* consequence of the design, not evidence against the score.

## 4. Gene sets (fixed before any model is fit)

**Primary repression-target set** (disjoint from zone markers; scored inversely —
higher target expression implies *lower* HDAC4 repressive activity):

    IBSP, SPP1, VEGFA, PANX3, SP7

**Secondary repression-target set** (full set as specified in the brief; reported as a
sensitivity analysis and flagged as partially circular with respect to zone assignment):

    COL10A1, MMP13, IBSP, SPP1, VEGFA, ALPL, PANX3, SP7

**Pathway-component set** (expression reported directly; **not** folded into the score):

    HDAC4, HDAC5, HDAC7, SIK1, SIK2, SIK3, CAMK4, CAMK2D, PPP2CA, PPP2R1A,
    PTH1R, PTHLH, MEF2C, RUNX2, SMAD2, SMAD3, HIF1A, EPAS1

**Housekeeping negative-control set** (scored by the identical procedure; must show no
zonal gradient):

    ACTB, GAPDH, TBP, RPL13A, PPIA, SDHA, B2M, HPRT1, PGK1, YWHAZ, RPLP0, UBC

**Scoring methods.** Two independent methods, both reported:

1. `decoupler` with a hand-curated MEF2C/RUNX2 regulon restricted to the target set
   (all targets weight +1, since HDAC4 relieves repression of MEF2C/RUNX2 targets).
2. Scaled mean expression via `scanpy.tl.score_genes` over the same set.

If the two methods disagree in direction, the disagreement is reported; neither is
picked as the winner.

**Sign convention, fixed now.** Both methods score *target expression*. HDAC4 repressive
activity is the negative of that. Reported "HDAC4 activity" = −(target score). The
zonal prediction is therefore: target score **increases** resting → hypertrophic;
equivalently inferred HDAC4 activity **decreases** resting → hypertrophic.

## 5. Primary comparison and the adequacy gate

**Primary comparison.** HDAC4 activity score in **resting-zone cells only**, younger
donors versus older donors, pseudobulk per donor.

**Adequacy gate (pre-specified, evaluated before the primary model is fit).** The
temporal analysis proceeds only if **both** hold:

- **(a)** Per-donor developmental age (or Tanner stage) is recorded, per donor, in the
  dataset metadata or an accompanying public record — not merely as a group-level range.
- **(b)** The recorded ages comprise ≥3 distinct values spanning ≥4 years, or ≥2 distinct
  Tanner stages.

If **(a)** fails, the primary comparison is **UNDEFINED** — not merely underpowered —
because the regressor does not exist. In that case no age model is fitted, no age
surrogate is invented, and donor index/accession order is **not** used as a proxy for
age. The finding is reported as the primary result.

If (a) holds but (b) fails, the analysis is reported as **UNINFORMATIVE** with the
effect size and interval shown and the inadequacy stated first.

**Zonal sanity gate.** The primary score must increase monotonically from resting to
hypertrophic (Spearman ρ > 0 across the four ordered zones, and hypertrophic > resting).
If this fails, the score is not working, nothing downstream is interpretable, and the
analysis stops and reports that.

## 6. Outcome definitions

| Verdict | Criterion |
|---|---|
| **Supporting** | RZ target-expression score increases with age (slope > 0), 95% CI excludes 0, **and** the boundary crossing point shifts toward the resting end with age |
| **Contradicting** | RZ target-expression slope on age < 0 with 95% CI excluding 0, or boundary crossing shifts toward the hypertrophic end |
| **Uninformative** | 95% CI includes 0; or the adequacy gate fails; or the zonal sanity gate fails; or any analysed zone has <100 cells |

A non-significant result is reported as non-significant. It is **not** described as a
"trend", "a suggestion of", or "approaching significance".

## 7. Pre-specified negative controls and guards

1. **Housekeeping control.** Score the housekeeping set by the identical procedure. It
   must show no zonal gradient. A gradient here invalidates the primary score.
2. **Label shuffle.** Shuffle zone labels within each sample and re-run. The zonal
   gradient must vanish. 1000 permutations; report the permutation null and empirical p.
3. **Cell-count reporting.** Report per-sample cell counts for every zone. Flag any zone
   with <100 cells; exclude flagged zone×sample cells from pseudobulk primary inference.
4. **Leave-one-out.** If an age model is fitted, refit leaving out each donor in turn to
   check whether any effect is driven by a single donor.
5. **Technical covariates.** Verify the score is not tracking library depth (total counts,
   n_genes) or mitochondrial fraction — report Spearman correlations at pseudobulk level.

## 8. Secondary analyses (pre-specified)

- **Boundary migration.** For each sample, locate the position along the differentiation
  ordering at which the target score crosses the midpoint between its resting and
  hypertrophic values. Test whether that crossing shifts toward the resting end with age.
  Conditional on the adequacy gate for the age test; the crossing position itself is
  computed and reported per sample regardless.
- **GH arm.** GH-treated versus vehicle organ culture, paired within donor, restricted to
  resting/progenitor cells. Pre-specified confound: 24 h culture is reported to lose the
  GP1 quiescent cluster entirely, so the contrast is restricted to GP2 and this is stated.
  This is the only interventional contrast in the data.
- **Pathway-component expression** by zone, reported as a table and heatmap, descriptively.

## 9. Deviations

Any deviation from this document is recorded in `docs/DEVIATIONS.md` with its reason.
