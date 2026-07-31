# Preregistration — Abaloparatide & final long-bone length

**Written before inspecting analytic results (only the target question and a handful of
already-known anchor-paper titles were known at authoring time). Timestamp basis: 2026-07-31.**

## Primary question (fixed)

In a growing organism with **open conventional long-bone growth plates**, does **intermittent
systemic abaloparatide** increase, decrease, or leave unchanged **final mineralized long-bone
length / adult height** vs no abaloparatide?

The primary endpoint is **stable, normally mineralized long-bone length at skeletal maturity**,
measured after treatment washout. Surrogates (cartilage thickness, COL10 level, open physis,
condyle length, short-term velocity) are explicitly **not** substitutes.

## Primary estimand

ΔL_final = L(mature, abaloparatide) − L(mature, control), measured at skeletal maturity after
washout and normal mineralization/remodeling.

Report as a continuous effect (percent of control long-bone length; cm of adult height) **and** as:

- P(ΔL_final > 0)   — any positive effect
- P(ΔL_final ≥ 1%)  (animals) / P(adult-height gain ≥ 1 cm) (humans) — meaningful benefit
- P(ΔL_final < 0)   — harm

## Secondary estimands (analyzed and reported SEPARATELY, never substituted for primary)

Short-term growth rate; height velocity; growth-plate thickness; resting-zone cell number/clonal
diversity; resting self-renewal vs commitment; BA/CA (maturation tempo); time to fusion;
proliferative output; entry into pre/hypertrophy; hypertrophic axial height/volume/dry-mass; matrix
production; primary-spongiosa/vascular-invasion/mineralization; bone geometry/alignment/mechanics;
**mandibular-condyle growth (separate category, never pooled with conventional physes)**;
fracture-callus/cartilage-repair (separate); **safety/mineral metabolism (separate from efficacy
sign)**.

## Evidence tiers (analyzed separately; combined only via explicit transportability weights)

- **T1** abaloparatide, conventional physis, growing organism, final mature length. (direct)
- **T2** abaloparatide, conventional physis, short-term growth/histology, not final maturity.
- **T3** abaloparatide, secondary cartilage / mandibular condyle.
- **T4** abaloparatide, fracture callus / cartilage repair.
- **T5** other intermittent PTH1R agonists in normal juvenile conventional physes.
- **T6** intermittent PTH1R agonists in disease/replacement (incl. pediatric hypoPTH).
- **T7** chronic developmental human genetics & transgenic models (boundary conditions only).
- **T8** molecular/cellular mechanism.

## Analysis rules (fixed in advance)

1. **Cells are never biological replicates.** Any interval from bootstrapping cells of one
   mouse/library is a *technical* sampling interval, labeled as such. Biological p-values require
   ≥ independent animals/experiments.
2. Cluster/annotate omics **without target genes**; freeze labels before querying Pth1r/effectors.
   Ihh may be used for annotation but any Ihh-dependent result is flagged **partially circular**.
3. No naive pooled meta-analysis across tiers; use transportability weighting with **broad
   sensitivity**. The final *sign* must not depend on one undisclosed weighting scheme.
4. Simulations (PK/PD, growth-plate ODE) are **constraint/sensitivity tools that expose what follows
   from assumptions**, never "evidence" and never "measured biology."
5. Efficacy sign is inferred **independently** of safety; safety reported separately.
6. Reviews used only to find primary sources; evidentiary conclusions rest on primary
   research / regulatory documents / registries / original datasets.

## Required final classifications

**Strict scientific classification** (direct-evidence only):
Demonstrated positive / Demonstrated neutral / Demonstrated negative / **Not identified by existing
direct evidence**.

**Forced best estimate:** Most likely positive / approximately neutral / negative — with posterior
probability, range across priors & transportability assumptions, strongest evidence for and against,
and the single missing observation that would most change it.

Decision thresholds: likely positive if P(ΔL>0) ≥ 0.75; likely negative if P(ΔL<0) ≥ 0.75;
otherwise likely neutral/indeterminate (also declared if the 80% interval lies mostly within ±1%).
If reasonable priors flip the category → verdict = "model-dependent / not identified" + a forced
research bet.

## Pre-declared hypotheses to test (not to assume)

- **Inverted-U**: both PTH1R LoF (Blomstrand) and constitutive gain (Jansen) shorten bones; net
  effect of an *intermittent partial* perturbation is interior and sign-ambiguous a priori.
- **Velocity≠final-length**: a short-term velocity increase can be spent on faster maturation
  (earlier fusion) with no adult-height gain, or even a decrease.
- **Intermittency**: pulsatile PTH1R agonism is mechanistically distinct from the tonic activation
  that causes Jansen dwarfism — must not be equated.

**The sign is NOT pre-decided.**
