#!/usr/bin/env python
"""
Monte-Carlo over abaloparatide candidate effects -> distribution of ΔL_final (%),
across 4 prior families x 3 developmental stages. Vectorized integrator (N draws at once).

Output: P(ΔL>0), P(ΔL>=1%), P(ΔL<0), median, 80/95% CI, plus P(short-term velocity up),
P(reserve preserved), and the joint P(velocity up AND no final-length gain).

SIMULATION = assumption-exposing constraint tool, NOT evidence. The headline scientific finding is
whether the SIGN is robust across prior families or prior-dominated (=> not identified).
"""
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 20260731
rng = np.random.default_rng(SEED)
N = 6000
dt = 0.1

BASE = dict(s0=0.55,a0=0.06,dR=0.010,s1=0.30,e=0.22,amp=4.0,kp=0.28,kc=0.45,
            kh=0.40,kv=0.45,hsize=1.00,matrix=1.00,K0=1.0)

# 13 candidate effects -> model nodes (documented in MODEL_CARD).
THETAS = ["surv","quies","renew","export","prolif","phdelay","ihh","matrix","hsize","vasc","spongiosa","osteo","desens"]
# node application: surv->dR-, quies->a0-, renew->s0+, export->e+, prolif->kp+, phdelay->kh-,
#   ihh (feedback on phdelay recovery)-> modulates phdelay effective, matrix->matrix+, hsize->hsize+,
#   vasc->kv+, spongiosa->kv+/matrix+ (mineralized replacement), osteo->matrix+ (bone-side), desens-> decays da.

def sample_family(fam):
    """Return dict of arrays (length N) for each theta + duty cycle d. Signs reflect the family's stance."""
    z = lambda mu, sd: rng.normal(mu, sd, N)
    if fam == "skeptical":
        base = {t: z(0.0, 0.04) for t in THETAS}
        d = rng.uniform(0.02, 0.08, N)
    elif fam == "neutral":
        base = {t: z(0.0, 0.22) for t in THETAS}   # symmetric, no directional prior
        d = rng.uniform(0.03, 0.12, N)
    elif fam == "mech_positive":
        # reserve preservation + maintained/added output (growth-promoting stance)
        base = dict(surv=z(0.25,0.15), quies=z(0.15,0.12), renew=z(0.20,0.15), export=z(0.10,0.12),
                    prolif=z(0.25,0.15), phdelay=z(0.20,0.15), ihh=z(0.10,0.10), matrix=z(0.20,0.15),
                    hsize=z(0.20,0.15), vasc=z(0.05,0.10), spongiosa=z(0.10,0.10), osteo=z(0.10,0.10),
                    desens=z(0.05,0.05))
        d = rng.uniform(0.04, 0.12, N)
    elif fam == "mech_negative":
        # reserve depletion / accelerated maturation / earlier fusion (growth-limiting stance)
        base = dict(surv=z(-0.15,0.12), quies=z(-0.20,0.15), renew=z(-0.10,0.12), export=z(0.20,0.15),
                    prolif=z(0.20,0.15), phdelay=z(0.05,0.12), ihh=z(-0.10,0.10), matrix=z(0.0,0.12),
                    hsize=z(0.0,0.12), vasc=z(0.20,0.15), spongiosa=z(0.15,0.12), osteo=z(0.05,0.10),
                    desens=z(0.20,0.12))
        d = rng.uniform(0.04, 0.12, N)
    base["desens"] = np.abs(base["desens"])   # desens >=0
    return base, np.clip(d, 0.01, 0.2)

STAGES = {   # (tc, w, B0, R0_init, R1_init) — early=reserve-rich/clock-far; late=reserve-low/clock-near
    "early_post_SOC":  dict(tc=135, w=8, B0=8.5, init=(0.95,0.15)),
    "mature_plate":    dict(tc=100, w=7, B0=6.0, init=(0.90,0.15)),
    "late_closing":    dict(tc=58,  w=5, B0=2.6, init=(0.45,0.10)),
}

def integrate(theta, d, stage, T=300.0):
    """Vectorized over N draws. theta: dict of arrays; d: array. Returns Lfinal array + peak early velocity."""
    tc, w, B0, (R0i, R1i) = stage["tc"], stage["w"], stage["B0"], stage["init"]
    n = int(T/dt)
    ones = np.ones(N)
    R0=ones*R0i; R1=ones*R1i; P=ones*0.10; PH=ones*0.05; H=ones*0.03; L=np.zeros(N); B=ones*B0
    fused = np.zeros(N, bool)
    des = theta["desens"]
    # precompute static base rates
    b = BASE
    Lprev = L.copy(); peakvel = np.zeros(N); half = tc  # early velocity window ~ before midpoint
    for i in range(n):
        t = i*dt
        g = 1.0/(1.0+np.exp((t-tc)/w))
        da = d*np.exp(-des*t)
        reserve_gone = (B<=1e-6) & ((R0+R1)<0.02)
        fused = fused | (g<0.02) | reserve_gone
        act = (~fused).astype(float)
        renew_ok = ((B>0)&(~fused)).astype(float)
        s0 = b["s0"]*np.exp(theta["renew"]*da)
        a0 = b["a0"]*np.exp(-theta["quies"]*da)
        dR = b["dR"]*np.exp(-theta["surv"]*da)
        e  = b["e"]*np.exp(theta["export"]*da)
        kp = b["kp"]*np.exp(theta["prolif"]*da)
        # ihh feedback: partially offsets phdelay over time (negative feedback recovers hypertrophy)
        kh = b["kh"]*np.exp(-(theta["phdelay"]-0.5*theta["ihh"])*da)
        kv = b["kv"]*np.exp((theta["vasc"]+0.5*theta["spongiosa"])*da)
        hsize = b["hsize"]*np.exp(theta["hsize"]*da)
        matrix = b["matrix"]*np.exp((theta["matrix"]+0.4*theta["osteo"]+0.3*theta["spongiosa"])*da)
        amp = b["amp"]
        renew_flux = s0*g*R0*(1-R0/b["K0"])*renew_ok
        B  = np.maximum(B - dt*renew_flux, 0)
        R0 = np.maximum(R0 + dt*(renew_flux - a0*R0*act - dR*R0), 0)
        R1 = np.maximum(R1 + dt*(a0*R0*act + b["s1"]*g*R1*renew_ok - e*R1*act - dR*R1), 0)
        P  = np.maximum(P  + dt*((e*R1*amp)*act + kp*g*P*act - b["kc"]*P), 0)
        PH = np.maximum(PH + dt*(b["kc"]*P - kh*PH), 0)
        H  = np.maximum(H  + dt*(kh*PH - kv*H), 0)
        L  = L + dt*(kv*H*hsize*matrix)*(~fused)
        if t < half:
            v = (L-Lprev)/dt
            peakvel = np.maximum(peakvel, v)
        Lprev = L.copy()
    return L, peakvel

results = {}
rows = []
for fam in ["skeptical","neutral","mech_positive","mech_negative"]:
    theta, d = sample_family(fam)
    zero = {t: np.zeros(N) for t in THETAS}
    for sname, stage in STAGES.items():
        Lc, vc = integrate(zero, d, stage)
        Ld, vd = integrate(theta, d, stage)
        dL = (Ld - Lc)/np.clip(Lc,1e-9,None)*100.0
        vel_up = vd > vc*1.001
        reserve_pres = Ld > Lc            # proxy handled below more specifically
        row = dict(family=fam, stage=sname,
                   P_dL_gt0=float(np.mean(dL>0)),
                   P_dL_ge1pct=float(np.mean(dL>=1.0)),
                   P_dL_lt0=float(np.mean(dL<0)),
                   P_dL_le_neg1pct=float(np.mean(dL<=-1.0)),
                   median_dL_pct=float(np.median(dL)),
                   dL_p10=float(np.percentile(dL,10)), dL_p90=float(np.percentile(dL,90)),
                   dL_p2_5=float(np.percentile(dL,2.5)), dL_p97_5=float(np.percentile(dL,97.5)),
                   P_velocity_up=float(np.mean(vel_up)),
                   P_velUp_and_no_final_gain=float(np.mean(vel_up & (dL<1.0))))
        rows.append(row); results[(fam,sname)] = dL
summary = pd.DataFrame(rows)
summary.to_csv("results/models/growth_deltaL_by_family_stage.csv", index=False)
pd.set_option("display.width",240,"display.max_columns",30)
print(summary.round(3).to_string(index=False))

# pooled across families (equal prior weight) for mature plate — the reference stage
mature = np.concatenate([results[(f,"mature_plate")] for f in ["skeptical","neutral","mech_positive","mech_negative"]])
pooled = dict(stage="mature_plate_pooled_equal_family_wt",
              P_dL_gt0=float(np.mean(mature>0)), P_dL_ge1pct=float(np.mean(mature>=1)),
              P_dL_lt0=float(np.mean(mature<0)), median=float(np.median(mature)),
              p2_5=float(np.percentile(mature,2.5)), p97_5=float(np.percentile(mature,97.5)))
json.dump(pooled, open("results/models/pooled_mature.json","w"), indent=2)
print("\nPOOLED (equal family weights, mature plate):")
for k,v in pooled.items(): print(f"  {k}: {v}")

# ---- figures ----
fig, axes = plt.subplots(1,3, figsize=(15,4.4), sharex=True)
for ax,(sname,_) in zip(axes, STAGES.items()):
    for fam,c in zip(["skeptical","neutral","mech_positive","mech_negative"],
                     ["#888888","#4C78A8","#54a24b","#e45756"]):
        ax.hist(np.clip(results[(fam,sname)],-8,8), bins=60, histtype="step", lw=1.6, color=c, label=fam, density=True)
    ax.axvline(0, color="k", lw=0.8); ax.axvline(1, color="orange", ls="--", lw=1, )
    ax.set_title(sname); ax.set_xlabel("ΔL_final (% of control)")
axes[0].set_ylabel("density"); axes[0].legend(fontsize=7)
fig.suptitle("Posterior ΔL_final by prior family & developmental stage (SIMULATION, not evidence)", fontsize=12)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig("results/figures/deltaL_posteriors.png", dpi=130)
print("\nSaved results/figures/deltaL_posteriors.png")
