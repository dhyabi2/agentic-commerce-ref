# Abaloparatide & final conventional long-bone length — evidence synthesis

**Best current answer: intermittent systemic abaloparatide is _not identifiable_ from direct evidence
for final conventional long-bone growth; forced best estimate = _approximately neutral_ (with a mild,
uncertain negative skew at the meaningful-effect threshold), at _low_ confidence.**

**Forced estimate: P(any benefit)=~0.36, P(meaningful benefit ≥1% animal / ≥1 cm human)=~0.13,
P(harm)=~0.64 any / ~0.30 meaningful (≤−1%); P(approximately neutral, |ΔL|<1%)=~0.57; median ΔL≈−0.3%,
95% interval [−10.6%, +5.3%]; range of medians across reasonable priors/stages = −1.7% (mechanism-
positive, early plate) to +0.44% (mechanism-positive, late-closing plate). The sign is
prior-dependent → formally model-dependent / not identified.**

---

## 1. Direct evidence (Tier 1–2): none

There is **no study** in which systemic abaloparatide is given to a normal growing organism with
**conventional** long-bone physes and final (or even short-term) long-bone length is measured. An
exhaustive Europe PMC sweep (~400 abaloparatide records) returns only: mandibular **condyle**
(secondary cartilage — PMID 36576115 systemic; PMID 37727017 **local** injection, +1.34 mm condylar
length) and **adult** metaphyseal/OVX bone (e.g. PMID 40254128, femur length not measured). This is a
genuine data gap, not an access failure.

→ **Strict scientific classification: NOT IDENTIFIED by existing direct evidence.**

## 2. Indirect evidence, by tier (never pooled naively; see EVIDENCE_LEDGER.csv)

- **T5 — closest analogue (Ogawa 2002, PMID 11862529):** hPTH(1-34) 80 µg/kg/d, 5×/wk × 3 wk in
  6-wk rats raised proximal-tibia growth-plate thickness, chondrocyte number and **longitudinal
  growth *rate*** — but only **short-term dynamic velocity**, no final length, no washout, and the
  effect was **absent in 15-wk adults**. Direction +, magnitude not retrieved. Velocity, not height.
- **T6 — pediatric hypoparathyroidism (Winer RCT PMID 20392870; cohorts PMC6298875, PMID 39883563):**
  PTH(1-34) **replacement** → normal height velocity, most reach adult height; randomized vs
  calcitriol shows **no growth difference**. This is **replacement-to-normal, not augmentation** in a
  normal child; serial bone ages are unavailable. It bounds PTH(1-34) as growth-**neutral at
  replacement dose**, and says nothing about supra-physiologic agonism.
- **T7 — human genetics (HUMAN_ALLELIC_SERIES.csv):** an **inverted-U**. Gs-arm loss (PTHLH-LoF/
  brachydactyly E2; Blomstrand; acrodysostosis PRKAR1A/PDE4D; PHP1a) → **advanced hypertrophy +
  premature fusion → short**. Constitutive receptor gain (Jansen) → **delayed, disorganized → also
  short**. Both chronic extremes shorten bones by opposite mechanisms; the naive "more signaling =
  more growth" is false. These are **boundary conditions**, not estimates of intermittent dosing.
- **T8 — mechanism:** PTH1R-Gs maintains resting-chondrocyte quiescence/survival (Chagin/Kronenberg
  2014, PMID 24781502); PTHrP+ resting cells are a **finite→self-renewing** stem pool that sets growth
  duration (Mizuhashi 2018 PMID 30401834; Newton 2019 PMID 30814736); estrogen **depletes** the
  reserve and **advances fusion** (Nilsson 2014, PMID 24708243). This is the structural basis for the
  reserve-limited model.
- **Comparator calibration:** GH raises attained adult height (+3.7 cm, PMID 15240584) — velocity
  *can* be banked. But **aromatase inhibitors raise *predicted* height +5.9 cm yet did NOT raise
  *attained* adult height** (PMID 16189252 / 31024444), and vosoritide raises velocity with final
  height still unknown (PMID 32891212). **Velocity/predicted gains frequently do not become final
  height** — the central caution for any PTH1R-agonist growth claim.

## 3. Receptor topology & reserve (Model 1) — audited from GSE288529 raw

Prior single-library audit (one 4-wk female mouse; descriptive only, cells ≠ biological replicates):
resting *Pth1r* pseudobulk **1.15** vs prehypertrophic **17.0 CP10K** (~15×), yet **63% of resting
cells detect Pth1r**. The operational-agonism sensitivity model gives median resting response ≈ **28%**
of prehypertrophic; P(resting ≥25/50/80% of prehyp) = **0.54 / 0.31 / 0.12**. **Conclusion:** low
resting *Pth1r* does **not** make direct resting-zone action implausible (resting is not silent), but
spare receptors **rarely fully erase** the gap. (`results/figures/pkpd_reserve_sensitivity.png`.)

## 4. Final-length model (Model 2) — the crux

A reserve-limited endochondral model (finite stem budget + developmental fusion clock) reproduces all
four boundary constraints (LoF short; tonic/Jansen short; estrogen-deficiency taller-but-bounded;
intermittent agonist raises early velocity) with one parameter set. Sweeping 13 candidate
abaloparatide effects (signs not assumed) across 4 prior families × 3 stages (6000 draws each;
`results/models/growth_deltaL_by_family_stage.csv`, `results/figures/deltaL_posteriors.png`):

- **skeptical / neutral:** ΔL centered ≈ 0 (median ~0; P(ΔL>0) ≈ 0.50); genuinely sign-ambiguous.
- **mechanism-positive:** paradoxically **negative** at early/mature stages (median −1.4 to −1.7%,
  P(<0) ≈ 0.77–0.81) — boosting output/flux spends the finite reserve faster → **earlier fusion**, so
  transient anabolic activity is **not banked** as length; turns **positive only at late-closing**
  (median +0.44%, P(>0) = 0.62, P(≥1%) = 0.38), when little reserve remains to deplete.
- **mechanism-negative:** mildly negative everywhere (median −0.4 to −0.5%).

The emergent finding — that even a growth-promoting prior can shorten final length via reserve
depletion — is exactly the aromatase-inhibitor lesson expressed mechanistically, and it is why
**velocity ≠ final length** dominates the answer.

## 5. Synthesis and forced estimate

No direct likelihood exists to anchor a Bayesian pooling (Tier 1 empty), so the posterior sign is
**prior-dominated**: it flips from mildly negative (mechanism-positive early/mature; mechanism-
negative) through ≈0 (neutral/skeptical) to mildly positive (mechanism-positive late-closing).
Under equal family weights at the mature plate: **P(any benefit)=0.36, P(≥1%)=0.13, P(any harm)=0.64
(meaningful ≥1%: 0.30), P(|ΔL|<1%)=0.57, median −0.28%, 95% CI [−10.6%, +5.3%].** Every accessible
indirect signal (Ogawa short-term-only; human replacement neutral; aromatase calibration; inverted-U
genetics) is concordant with **≈ neutral final length**, and leave-one-tier removal of any indirect
tier does not change the category.

### Required classifications
- **Strict:** **Not identified** by existing direct evidence.
- **Forced best estimate:** **Most likely approximately neutral** for final conventional long-bone
  length, with a mild and uncertain skew toward net shortening at the meaningful-effect threshold
  (P(meaningful harm 0.30) > P(meaningful benefit 0.13)). Because reasonable priors flip the sign, the
  formal verdict is **model-dependent / not identified**.
- **Forced research bet (as required):** if compelled to bet on the outcome of the decisive experiment
  in a mid-productive plate, I bet **approximately neutral, most-likely-slightly-negative** on final
  length despite a **likely transient velocity increase** — because in a reserve-limited system
  transient stimulation tends to be spent on tempo, not banked as height, and abaloparatide's
  RG-biased transient signaling argues against a large sustained anabolic drive. The **strongest
  positive scenario** is a **late-pubertal adjunct** (spend residual reserve while the clock is
  closing), and the **strongest negative scenario** is early/mid dosing that advances maturation and
  fuses earlier.

**Strongest evidence FOR a benefit:** Ogawa short-term velocity + PTH1R-Gs maintains resting
quiescence (could preserve reserve) + intermittency avoids tonic Jansen-like drag.
**Strongest evidence AGAINST:** inverted-U genetics (both extremes short) + aromatase-inhibitor
calibration (velocity/predicted gains not banked) + reserve-depletion dynamics + human replacement
data showing normalization, not augmentation.
**Single most decisive missing observation:** final mature long-bone length (+ serial bone age /
time-to-fusion) in normal juveniles given intermittent systemic abaloparatide with washout to
maturity (see DATA_REQUEST.md §4).

## 6. Safety (reported separately from efficacy sign)
Adult-rat carcinogenicity (FDA label §13.1): dose-dependent osteosarcoma (0–2% → 87% M/62% F at
28× human exposure) — a **bone-proliferation/tumor** signal in mature rats, **not** a physeal or
long-bone-length finding. Whether any **juvenile** toxicology or physeal histology was submitted is
unconfirmed (DATA_GAPS §3). This safety signal did **not** enter the efficacy inference.

## 7. Limitations
Single-library omics (no biological replicate); direct tier empty; models are assumption-exposing
simulations, not measured biology; several regulatory PDFs and paywalled kinetic constants not
retrieved; condyle and chronic-genetics evidence deliberately down-weighted per preregistration.

---
*See `docs/PREREGISTRATION.md`, `docs/MODEL_CARD.md`, `docs/DATA_GAPS.md`, `docs/DATA_REQUEST.md`,
`results/tables/*.csv`, `results/figures/*.png`. Reproduce with `bash run_all.sh`.*
