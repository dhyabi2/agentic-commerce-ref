# Abaloparatide & final long-bone length — definitive evidence synthesis

**Question:** In a growing organism with open **conventional** long-bone physes, does intermittent
systemic **abaloparatide** change **final mineralized long-bone length / adult height**?

**Answer:** *Not identifiable from direct evidence; forced best estimate ≈ **approximately neutral**
(mild, uncertain negative skew), low confidence.* See `REPORT.md` (full) and `EXECUTIVE_SUMMARY.md`.

Forced estimate: P(any benefit)=0.36 · P(meaningful benefit ≥1%/≥1 cm)=0.13 · P(harm)=0.64 any /
0.30 meaningful · median ΔL≈−0.3% · 95% CI [−10.6%, +5.3%]. Sign is prior-dependent → not identified.

## Reproduce
```bash
bash run_all.sh          # builds tables, runs models, writes figures (~40 s; no network needed)
```
Environment pins + OS + seed + input checksum in `environment.txt` (seed 20260731).

## Layout
```
docs/PREREGISTRATION.md   estimands & rules fixed before analysis
docs/MODEL_CARD.md        equations, priors, calibration, non-identifiability, validation
docs/DATA_GAPS.md         what is missing, ranked by impact
docs/DATA_REQUEST.md      exact data to request (NIH bone ages, FDA/EMA PDFs, atlas raw, the experiment)
src/                      builders + 3 models (evidence ledger, PK/PD reserve, growth-plate + MC, figures)
results/tables/           EVIDENCE_LEDGER.csv, HUMAN_ALLELIC_SERIES.csv, CLAIMS_AUDIT.csv
results/models/           growth_deltaL_by_family_stage.csv, pooled_mature.json
results/sensitivity/      pkpd reserve summaries
results/figures/          evidence_map, pkpd_reserve_sensitivity, deltaL_posteriors
REPORT.md / EXECUTIVE_SUMMARY.md
```

## Method in one paragraph
Evidence was harvested from primary literature, registries and regulatory documents (verified PMIDs/
URLs in `EVIDENCE_LEDGER.csv`), stratified into 8 non-pooled tiers by directness. **Tier 1–2 (direct
abaloparatide, conventional physis, final length) are empty** — a genuine gap. Three assumption-
exposing models (labeled *not evidence*) test (1) whether low resting-zone *Pth1r* precludes direct
action, (2) how 13 candidate drug effects map to final length in a **reserve-limited** growth plate,
and (3) whether the sign is robust across prior families. It is not — hence *not identified*, with a
forced approximately-neutral best estimate. Efficacy and safety are analyzed separately; cells are
never treated as biological replicates.

*Scientific research synthesis — not a treatment or dosing recommendation.*
