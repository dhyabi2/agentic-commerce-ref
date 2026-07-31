# MODEL CARD — abaloparatide / final long-bone length

Three models. **All are assumption-exposing constraint/sensitivity tools, not evidence and not
measured biology.** Their job is to show *what follows from* the stated priors and the verified
pharmacology/topology — especially whether the sign of the final-length effect is robust or
prior-dominated.

Seeds fixed (`SEED = 20260731`). Deterministic given seed. Python 3.11, numpy/scipy/pandas/matplotlib
(pins in `environment.txt`).

---

## Model 1 — PK/PD receptor reserve (`src/pkpd_receptor_reserve.py`)

**Purpose:** can spare receptors erase the ~15× resting-vs-prehypertrophic *Pth1r* mRNA gap
(GSE288529), making direct resting-zone action plausible?

**A. Operational (Black–Leff) agonism.** Response fraction per zone
`E = Emax·(τ·occ)^n / (1+(τ·occ)^n)`, with transduction ratio τ ∝ receptor protein density.
Receptor mRNA per zone → protein via a shared sublinear exponent `p∈U(0.5,1)` + lognormal zone noise;
τ anchored to the prehypertrophic (highest) zone, `τ_ref` log-uniform 0.3–100 (reserve); cartilage
occupancy `occ∈(0,1)` skewable low; slope `n∈U(0.8,1.5)`. 200 000 draws.

**B. Compact cAMP ODE.** `d[cAMP]/dt = kcat·R_active·occ·pulse(t) − kpde·cAMP`; intermittent pulse
2–20 min (abaloparatide RG-selective → short); `kpde` log-uniform (t½ 0.7–14 min). Reports cAMP-AUC
and time-above-threshold per zone.

**Key priors are deliberately WIDE** because cartilage PK and zonal receptor protein are unmeasured.

**Headline outputs:** median resting/prehyp response ≈ **0.28**; P(resting ≥25%/≥50%/≥80% of prehyp)
= **0.54 / 0.31 / 0.12**; resting effectively silent (<5% Emax) in **0.38**. Gap erased (≥80%) only
under high reserve **and** high occupancy (0.52 vs 0.05 otherwise). Interpretation: resting is **not
signaling-silent**, but spare receptors **rarely fully erase** the gap.

---

## Model 2 — Dynamic growth-plate model (`src/growth_plate_model.py`, MC in `src/run_growth_synthesis.py`)

**States:** R0 (root resting) · R1 (active PTHrP+ resting) · P (transit-amplifying) · PH
(prehypertrophic) · H (hypertrophic) · L (cumulative mineralized length = primary endpoint).

**Two binding constraints (the crux of velocity ≠ final length):**
- **CLOCK:** developmental senescence `sen(t)=1/(1+exp((t−tc)/w))` closes the window (fusion) near tc
  (estrogen-driven; PMID 24708243).
- **RESERVE:** a **finite** stem-renewal budget B0. R0 self-renewal consumes B0; when spent, R0
  cannot renew and the pool drains. Faster flux spends reserve sooner → earlier stop.

`dL/dt = kv·H·hsize·matrix` while unfused; **elongation stops at fusion**. Final length =
`min(reserve-limited, clock-limited)` accrual, so an intervention lengthens the bone only if it
relaxes the actually-binding constraint. Transit pool kept input-driven (kp<kc) so total output is
reserve-throughput-limited, not self-exploding.

**Drug = intermittent PTH1R agonism**: 13 candidate effects (signs NOT assumed) applied as log-fold
rate multipliers, gated by duty cycle `d` (0.01–0.2, hours/day — the intermittency knob) and optional
desensitization. Mapping of the 13 → model nodes:

| candidate effect | node |
|---|---|
| 1 resting survival | dR− |
| 2 resting quiescence | a0− |
| 3 self-renewal/commitment | s0+ |
| 4 export into columns | e+ |
| 5 proliferative division | kp+ |
| 6 prehypertrophic delay | kh− |
| 7 IHH suppression/recovery | partial offset of kh (feedback) |
| 8 matrix synthesis | matrix+ |
| 9 hypertrophic size | hsize+ |
| 10 vascular invasion | kv+ |
| 11 primary spongiosa/mineralization | kv+ / matrix+ |
| 12 osteoblast bone formation | matrix+ |
| 13 receptor downregulation/desensitization | decays drug activity over t |

**Prior families:** skeptical (≈0 effects); neutral (symmetric, no directional prior);
mechanism-positive (reserve-preservation + added output); mechanism-negative (reserve depletion /
accelerated maturation). **Stages:** early-post-SOC (reserve-rich, clock-far); mature-plate;
late-closing (reserve-low, clock-near). N=6000 draws each; control = all effects 0.

**Validation (`validate()` — all pass with one consistent parameter set):**
- PTH1R LoF (↑kh, ↓resting maintenance) → **shorter** (premature reserve exhaustion). ✓
- Tonic constitutive gain / Jansen (strong phdelay, d=1) → **shorter** (cells stuck pre-hypertrophic,
  fail to mineralize before fusion). ✓
- Estrogen deficiency (clock delayed, tc↑) → **taller but bounded** by finite reserve. ✓
- Intermittent agonist → **raised early velocity**. ✓ (reproduces Ogawa 2002 direction)

**Non-identifiability (stated plainly):** the ΔL_final **sign is prior-family- and stage-dependent**
(mechanism-positive is negative at early/mature but positive at late-closing; neutral/skeptical ≈ 0).
No direct data pin the effect parameters, so the posterior sign is **not identified** by data; the
pooled estimate is a construction under equal family weights, reported as such.

---

## Model 3 — Evidence synthesis (transportability, in `REPORT.md` §synthesis)

Because Tier-1/2 (direct) evidence is **empty**, no same-drug final-length likelihood exists to anchor
a Bayesian pooling. The synthesis is therefore explicit: direct status = *not identified*; the forced
estimate is the Model-2 posterior modulated by weak, concordant indirect signals (Ogawa short-term +;
human replacement neutral final height; aromatase-inhibitor calibration that velocity/predicted gains
often are not banked; inverted-U genetics). We run the sign across prior families rather than pretend
one weighting is objective. Leave-one-tier: removing any single indirect tier does not move the
category (all point to ≈ neutral final length); removing the (empty) direct tier is what would move it.

## What would most reduce uncertainty (value of information)
**A single direct measurement:** final mature long-bone length (+ serial bone age / time-to-fusion)
in normal juvenile animals given intermittent systemic abaloparatide, with washout to skeletal
maturity. This alone would populate the empty Tier 1 and collapse most of the posterior width.
