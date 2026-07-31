# Model card — seed digital twin (SIMULATION, not evidence)

The warehouse ships a **seed** digital twin reused from the companion abaloparatide analysis. It is an
assumption-exposing constraint tool. **It is not evidence and not measured biology**; outputs are what
follows from priors.

## Models present
1. **Reserve-limited growth-plate ODE** (`src/modeling/growth_plate_ode.py`,
   `src/modeling/twin_intervention_sweep.py`). States R0/R1/P/PH/H → cumulative length L. Two binding
   constraints: a **finite stem-renewal budget** and a **developmental fusion clock**; final length =
   min(reserve-limited, clock-limited). A generic intervention is a 13-effect log-fold vector applied
   at an intermittent duty cycle. Reproduces four boundary constraints (LoF→short; tonic gain→short;
   estrogen-deficiency→taller-but-bounded; intermittent agonist→↑velocity). Sweep over 4 prior families
   × 3 stages → ΔL_final distribution (`results/models/growth_deltaL_by_family_stage.csv`,
   `results/figures/deltaL_posteriors.png`). Key emergent result: **transient velocity gains are often
   not banked as final length** (reserve depletion → earlier fusion) — the mechanistic form of the
   surrogate-calibration lesson.
2. **PK/PD receptor-reserve** (`src/pkpd/pkpd_receptor_reserve.py`). Operational (Black–Leff) agonism +
   compact cAMP ODE; tests whether low resting-zone receptor mRNA precludes response
   (`results/figures/pkpd_reserve_sensitivity.png`).

## Not implemented (specified in §16–17 of the brief; see LIMITATIONS)
Spatial reaction-diffusion; agent-based lineage; FEniCS finite-element mechanical twin; full Bayesian
(PyMC/Stan) calibration with Sobol indices. The ODE sweep uses transparent Monte-Carlo over prior
families instead of PyMC.

## Calibration controls available in the corpus (for future full calibration)
estrogen/aromatase deficiency (CYP19A1/ESR1 → tall, unfused), FGFR3 gain/loss, NPR2 gain/loss,
GH deficiency/replacement, PTH1R/PTHLH series, glucocorticoid, known final-height trials — all in
`genes`/`surrogate_calibration`/`studies`.
