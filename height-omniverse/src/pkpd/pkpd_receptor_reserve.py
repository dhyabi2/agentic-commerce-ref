#!/usr/bin/env python
"""
PK/PD receptor-reserve + cAMP-dynamics sensitivity model.

Question addressed (Section 5): across parameter sets consistent with known pharmacology and the
GSE288529 receptor-mRNA topology, is low resting-zone Pth1r sufficient to make DIRECT resting-zone
action implausible, or can spare receptors plausibly erase the transcript gap?

This is a CONSTRAINT / SENSITIVITY tool. It is NOT measured biology. Outputs are "what follows from
assumptions", with wide transparent priors because cartilage PK and zonal receptor protein are not
directly measured.

Two coupled analyses:
  A. Black-Leff operational agonism -> steady response fraction per zone, reserve captured by tau.
  B. Compact cAMP ODE (production ~ active receptor; PDE decay) -> cAMP AUC & time-above-threshold
     per zone for an intermittent (pulsatile) ligand exposure. Abaloparatide is modeled as
     RG-selective -> shorter effective signaling pulse than PTH(1-34) (sensitivity-tested).

Receptor mRNA input = GSE288529 per-state Pth1r pseudobulk (CP10K), ONE library/ONE animal ->
descriptive relative abundances only; no biological CIs generated here.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20260731
rng = np.random.default_rng(SEED)
N = 200_000

# --- Receptor mRNA per zone (CP10K), from prior GSE288529 audit (single library) ---
mrna = {"Resting": 1.154, "Proliferative": 2.797, "Prehypertrophic": 16.999, "Hypertrophic": 16.150}
zones = list(mrna)
ref = "Prehypertrophic"                      # reference zone for tau anchoring (highest receptor)

# =====================================================================================
# Priors (wide, transparent). log-uniform unless noted.
# =====================================================================================
def loguni(lo, hi, n=N):
    return 10 ** rng.uniform(np.log10(lo), np.log10(hi), n)

# mRNA -> protein: shared sublinear exponent p (0.5..1.0) allows the gap to compress at protein level;
# plus independent lognormal zone noise (protein/mRNA ratio varies by zone, unmeasured).
p_exp   = rng.uniform(0.5, 1.0, N)
# transduction ratio at the reference (prehyp) zone: reserve. tau_ref ~ 0.3 (poorly coupled) .. 100 (huge reserve)
tau_ref = loguni(0.3, 100.0)
# ligand occupancy (fraction) in cartilage, same across zones (well-mixed within thin plate).
# Driven by unknown cartilage [A]/KA -> occupancy spans deep-linear to near-saturating.
occ     = rng.uniform(0.0, 1.0, N) ** rng.uniform(0.5, 3.0, N)   # skewable toward low occupancy
occ     = np.clip(occ, 1e-4, 0.999)
# operational slope
nH      = rng.uniform(0.8, 1.5, N)
Emax    = 1.0

def zone_tau(zone):
    # protein density relative to reference, with sublinear map + lognormal zone noise
    rel = (mrna[zone] / mrna[ref]) ** p_exp
    znoise = np.exp(rng.normal(0.0, 0.5, N))          # ~1.6x geometric SD
    return tau_ref * rel * znoise

def op_response(tau):
    x = (tau * occ) ** nH
    return Emax * x / (1.0 + x)

resp = {z: op_response(zone_tau(z)) for z in zones}
ratio_rest_pre = resp["Resting"] / np.clip(resp["Prehypertrophic"], 1e-9, None)

# =====================================================================================
# A. Operational-agonism results
# =====================================================================================
qs = [2.5, 25, 50, 75, 97.5]
rowsA = []
for z in zones:
    r = resp[z]
    rowsA.append({"zone": z, "mrna_CP10K": mrna[z],
                  **{f"resp_p{q}": np.percentile(r, q) for q in qs}})
opA = pd.DataFrame(rowsA).set_index("zone")

P25 = np.mean(ratio_rest_pre >= 0.25)
P50 = np.mean(ratio_rest_pre >= 0.50)
P80 = np.mean(ratio_rest_pre >= 0.80)
Pnull_resting = np.mean(resp["Resting"] < 0.05)   # resting response essentially silent

# Regime attribution: is the gap erased mainly when reserve high AND occupancy high?
high_reserve = tau_ref > 10
high_occ = occ > 0.5
erased = ratio_rest_pre >= 0.80

summaryA = {
    "P(resting >= 25% of prehyp response)": float(P25),
    "P(resting >= 50% of prehyp response)": float(P50),
    "P(resting >= 80% of prehyp response)": float(P80),
    "P(resting response effectively silent <5% Emax)": float(Pnull_resting),
    "median resting/prehyp response ratio": float(np.median(ratio_rest_pre)),
    "P(gap erased | high reserve & high occupancy)": float(np.mean(erased[high_reserve & high_occ])),
    "P(gap erased | low reserve OR low occupancy)": float(np.mean(erased[~(high_reserve & high_occ)])),
}

# =====================================================================================
# B. Compact cAMP ODE (per zone), intermittent pulse.
#   d[cAMP]/dt = kcat * Ractive(t) - (kpde) * cAMP
#   Ractive(t) = Rdens_zone * occ * pulse(t), pulse = 1 during ligand window, else 0.
#   Abaloparatide RG-selective -> shorter effective signaling window than PTH(1-34).
#   We report cAMP AUC per zone and time above a response threshold, relative to prehyp.
# =====================================================================================
def camp_auc(Rdens_rel, pulse_min, kpde_permin, tmax=240.0, dt=0.5):
    t = np.arange(0, tmax, dt)
    camp = np.zeros_like(Rdens_rel)  # vectorized over draws; start 0
    auc = np.zeros_like(Rdens_rel)
    above = np.zeros_like(Rdens_rel)
    thr = 0.2  # relative threshold
    prod0 = Rdens_rel * occ
    for ti in t:
        active = np.where(ti < pulse_min, prod0, 0.0)
        camp = camp + dt * (active - kpde_permin * camp)
        auc += camp * dt
        above += (camp > thr) * dt
    return auc, above

# pulse length prior: abaloparatide (RG) shorter; test 2-20 min effective signaling window
pulse_min = rng.uniform(2.0, 20.0, N)
kpde = loguni(0.05, 1.0)     # cAMP decay rate (per min); t1/2 ~ 0.7..14 min
Rrel = {z: (mrna[z]/mrna[ref])**p_exp * np.exp(rng.normal(0,0.5,N)) for z in zones}
aucB, aboveB = {}, {}
for z in zones:
    a, ab = camp_auc(Rrel[z], pulse_min, kpde)
    aucB[z] = a; aboveB[z] = ab
auc_ratio_rest_pre = aucB["Resting"] / np.clip(aucB["Prehypertrophic"], 1e-9, None)

summaryB = {
    "median cAMP-AUC resting/prehyp": float(np.median(auc_ratio_rest_pre)),
    "P(cAMP-AUC resting >= 50% prehyp)": float(np.mean(auc_ratio_rest_pre >= 0.5)),
    "P(cAMP-AUC resting >= 80% prehyp)": float(np.mean(auc_ratio_rest_pre >= 0.8)),
}

# =====================================================================================
# Save
# =====================================================================================
opA.to_csv("results/sensitivity/pkpd_operational_response_by_zone.csv")
allsum = {**{"[A]operational " + k: v for k, v in summaryA.items()},
          **{"[B]cAMP " + k: v for k, v in summaryB.items()}}
pd.Series(allsum).to_csv("results/sensitivity/pkpd_reserve_summary.csv", header=["value"])

print("=== A. Operational response by zone (percentiles) ===")
print(opA.round(3).to_string())
print("\n=== A. Reserve/occupancy erasure of the resting-prehyp gap ===")
for k, v in summaryA.items():
    print(f"  {k:55s}: {v:.3f}")
print("\n=== B. cAMP-dynamics (intermittent pulse) ===")
for k, v in summaryB.items():
    print(f"  {k:40s}: {v:.3f}")

# --- Figure: sensitivity surface P(resting>=80% prehyp) over reserve x occupancy ---
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
# heatmap of erasure probability
tb = np.log10(tau_ref); ob = occ
tg = np.linspace(tb.min(), tb.max(), 22); og = np.linspace(0, 1, 22)
H = np.full((len(og)-1, len(tg)-1), np.nan)
for i in range(len(og)-1):
    for j in range(len(tg)-1):
        m = (ob>=og[i])&(ob<og[i+1])&(tb>=tg[j])&(tb<tg[j+1])
        if m.sum()>50: H[i,j] = np.mean(erased[m])
im = axes[0].imshow(H, origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=1,
                    extent=[tg[0], tg[-1], 0, 1])
axes[0].set_xlabel("log10 receptor reserve tau (prehyp ref)"); axes[0].set_ylabel("cartilage occupancy")
axes[0].set_title("P(resting reaches >=80% of prehyp response)")
fig.colorbar(im, ax=axes[0], shrink=0.85)
# distribution of resting/prehyp response ratio
axes[1].hist(ratio_rest_pre, bins=60, color="#4C78A8", alpha=0.9)
for q,c,l in [(0.25,"#e45756","25%"),(0.5,"#54a24b","50%"),(0.8,"#f58518","80%")]:
    axes[1].axvline(q, color=c, ls="--", lw=1.5, label=f"{l} of prehyp")
axes[1].set_xlabel("resting / prehyp operational response ratio"); axes[1].set_ylabel("draws")
axes[1].set_title("Resting-zone response relative to prehypertrophic"); axes[1].legend(fontsize=8)
fig.suptitle("PK/PD receptor-reserve sensitivity (assumption-exposing, not measured biology)", fontsize=11)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("results/figures/pkpd_reserve_sensitivity.png", dpi=130)
print("\nSaved results/figures/pkpd_reserve_sensitivity.png")
