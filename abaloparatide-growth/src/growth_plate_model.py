#!/usr/bin/env python
"""
Dynamic growth-plate model: local signaling -> FINAL long-bone length (Section 6).

Compartments (reserve-limited endochondral column):
  R0 : deep/root resting progenitors (finite, senescence-gated self-renewal)
  R1 : active PTHrP+ resting progenitors
  P  : transit-amplifying proliferative cells
  PH : prehypertrophic cells
  H  : hypertrophic cells
  L  : cumulative mineralized length (the PRIMARY endpoint)

Core biology encoded:
 - Total elongation is RESERVE-LIMITED: R0 self-renewal is gated by a developmental senescence
   clock sen(t) (estrogen-driven). When the reserve exhausts, the plate fuses and growth stops.
 - dL/dt = (rate hypertrophic cells are replaced by bone) x (hypertrophic axial size) + matrix term.
 - This makes VELOCITY and FINAL LENGTH dissociable: boosting flux raises velocity but can spend
   the reserve faster (earlier fusion), so final length need not rise (calibrated to the aromatase-
   inhibitor lesson: predicted/velocity gains often do not become adult-height gains).

Drug = intermittent (pulsatile) PTH1R agonism, applied at duty cycle d (small: hours/day), with
optional desensitization. 13 candidate effects, signs NOT assumed; priors set per mechanism family.

THIS IS A SIMULATION: it exposes what follows from assumptions. It is not evidence/measured biology.
Calibrated qualitatively to boundary constraints (see validate()).
"""
import numpy as np
import pandas as pd

SEED = 20260731
rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------ baseline params
BASE = dict(
    s0=0.55,     # R0 self-renewal (gated by sen)
    a0=0.06,     # R0 -> R1 activation
    dR=0.010,    # resting death
    s1=0.30,     # R1 self-maintenance (gated)
    e=0.22,      # R1 -> P export
    amp=4.0,     # transit amplification factor on export
    kp=0.28,     # P proliferation (gated) -- kept < kc so P is input-driven, not self-exploding
    kc=0.45,     # P -> PH commitment (transit pool drains to its input)
    kh=0.40,     # PH -> H transition
    kv=0.45,     # H -> bone (vascular invasion/mineralization)
    hsize=1.00,  # hypertrophic axial size (length per replaced H)
    matrix=1.00, # matrix displacement multiplier on dL
    K0=1.0,      # R0 carrying capacity
)

def sen(t, tc, w):
    """developmental senescence clock: 1 (young) -> 0 (fused). tc=midpoint, w=width."""
    return 1.0 / (1.0 + np.exp((t - tc) / w))

def simulate(p, tc=100.0, w=7.0, T=260.0, dt=0.05, drug=None, B0=6.0):
    """Integrate to fusion. Final length = min(reserve-limited, clock-limited) accrual.

    Two binding constraints:
      * CLOCK: developmental senescence sen(t) closes the window near t=tc; elongation stops (fusion).
      * RESERVE: a FINITE stem-renewal budget B0. R0 self-renewal consumes budget; when exhausted,
        R0 can no longer renew and the pool drains. Faster flux spends reserve sooner (earlier stop).
    An intervention lengthens the bone only if it relaxes the ACTUALLY-binding constraint; transient
    velocity gains trade against duration.
    drug: dict of theta log-effects + duty cycle d + desens.
    """
    R0, R1, P, PH, H, L = 0.9, 0.15, 0.10, 0.05, 0.03, 0.0
    B = B0                                   # finite renewal budget (stem divisions remaining)
    n = int(T / dt)
    def eff(name):
        return drug.get(name, 0.0) if drug else 0.0
    d = drug.get("d", 0.0) if drug else 0.0
    des = drug.get("desens", 0.0) if drug else 0.0
    Ls = []; fused = False
    for i in range(n):
        t = i * dt
        g = sen(t, tc, w)                    # developmental clock gate
        da = d * np.exp(-des * t) if drug else 0.0
        # fusion when clock closes OR reserve+pool exhausted
        reserve_gone = (B <= 1e-6) and (R0 + R1) < 0.02
        if (g < 0.02) or reserve_gone:
            fused = True
        # --- drug-modulated rates (theta>0 raises the named process) ---
        s0 = p["s0"] * np.exp(eff("renew") * da)
        a0 = p["a0"] * np.exp(-eff("quies") * da)
        dR = p["dR"] * np.exp(-eff("surv") * da)
        e  = p["e"]  * np.exp(eff("export") * da)
        kp = p["kp"] * np.exp(eff("prolif") * da)
        kh = p["kh"] * np.exp(-eff("phdelay") * da)
        kv = p["kv"] * np.exp(eff("vasc") * da)
        hsize = p["hsize"] * np.exp(eff("hsize") * da)
        matrix = p["matrix"] * np.exp(eff("matrix") * da)
        amp = p["amp"] * np.exp(eff("ampout") * da)

        renew_active = (B > 0) and (not fused)
        prod = 0.0 if fused else 1.0         # production halts at fusion
        # self-renewal consumes finite budget
        renew_flux = s0 * g * R0 * (1 - R0/p["K0"]) * (1.0 if renew_active else 0.0)
        dB  = -renew_flux
        dR0 = renew_flux - a0*R0*prod - dR*R0
        dR1 = a0*R0*prod + p["s1"]*g*R1*(1.0 if renew_active else 0.0) - e*R1*prod - dR*R1
        dP  = (e*R1*amp)*prod + kp*g*P*prod - kc0(p)*P
        dPH = kc0(p)*P - kh*PH
        dH  = kh*PH - kv*H
        dL  = (kv*H * hsize * matrix) if not fused else 0.0   # elongation stops at fusion
        B=max(B+dt*dB,0); R0+=dt*dR0; R1+=dt*dR1; P+=dt*dP; PH+=dt*dPH; H+=dt*dH; L+=dt*dL
        R0=max(R0,0); R1=max(R1,0); P=max(P,0); PH=max(PH,0); H=max(H,0)
        Ls.append(L)
        if fused and (H < 1e-4) and (PH < 1e-4):
            break
    Ls = np.array(Ls)
    vel = np.gradient(Ls, dt) if len(Ls) > 2 else np.array([0.0])
    return dict(Lfinal=L, t=np.arange(len(Ls))*dt, L=Ls, vel=vel)

def kc0(p):
    return p["kc"]

# ------------------------------------------------------------------ validation
def validate():
    """Qualitative boundary checks. Returns dict of pass/fail + numbers."""
    base = simulate(BASE)
    Lb = base["Lfinal"]
    # PTH1R LoF: lose PTHrP delay -> faster hypertrophy (kh up), less resting maintenance -> shorter
    lof = BASE.copy(); lof["kh"] *= 1.8; lof["s1"] *= 0.6; lof["a0"] *= 1.5
    Llof = simulate(lof)["Lfinal"]
    # Constitutive/ tonic excess (Jansen): strong tonic PH->H block, d=1 (no intermittency)
    # -> cells stuck pre-hypertrophic, fail to mineralize before fusion -> disorganized SHORT
    jansen = simulate(BASE, drug=dict(phdelay=2.0, d=1.0, desens=0.0))["Lfinal"]
    # Estrogen deficiency: developmental clock delayed (later fusion) -> taller but BOUNDED by finite reserve
    estro = simulate(BASE, tc=140.0)["Lfinal"]
    # Short-term velocity bump from intermittent agonist (Ogawa): peak velocity rises early
    drug_v = simulate(BASE, drug=dict(prolif=0.35, phdelay=0.3, d=0.08, desens=0.0))
    base_pk = base["vel"][: len(base["vel"])//2].max()
    drug_pk = drug_v["vel"][: len(drug_v["vel"])//2].max()
    checks = {
        "baseline_Lfinal": Lb,
        "PTH1R_LoF_shorter_than_baseline": Llof < Lb, "LoF_L": Llof,
        "Jansen_tonic_shorter_than_baseline": jansen < Lb, "Jansen_L": jansen,
        "estrogen_deficiency_taller_(later_fusion)": estro > Lb, "estro_L": estro,
        "intermittent_agonist_raises_early_velocity": drug_pk > base_pk,
        "base_peakvel": base_pk, "drug_peakvel": drug_pk,
    }
    return checks

if __name__ == "__main__":
    import json
    v = validate()
    print("=== MODEL VALIDATION (qualitative boundary constraints) ===")
    for k, val in v.items():
        print(f"  {k:48s}: {val}")
    ok = all(v[k] for k in ["PTH1R_LoF_shorter_than_baseline","Jansen_tonic_shorter_than_baseline",
                            "estrogen_deficiency_taller_(later_fusion)","intermittent_agonist_raises_early_velocity"])
    print("ALL CORE CONSTRAINTS PASSED:", ok)
