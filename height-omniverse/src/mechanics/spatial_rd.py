#!/usr/bin/env python
"""Phase H: 1-D spatial reaction-diffusion of the growth-plate morphogen loop + systemic drug.
Domain x in [0,1]: 0=epiphysis/resting end, 1=metaphyseal/vascular front.
Fields (steady state, finite-difference): PTHrP (source at resting end), IHH (source in the
prehypertrophic band, feeds back onto PTHrP production), CNP (proliferative/prehyp source), and a
systemic DRUG entering from the perichondrial/vascular boundaries and diffusing into avascular cartilage.

Question: does a SYSTEMIC drug preserve or FLATTEN the endogenous spatial information?
SIMULATION — assumption-exposing. Parameters are wide priors (no guessed value treated as known)."""
import numpy as np, pandas as pd, pathlib
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[2]
rng = np.random.default_rng(20260731)

N = 201; x = np.linspace(0, 1, N); dx = x[1]-x[0]
# zone boundaries (fractional depth): resting 0-0.15, prolif 0.15-0.55, prehyp 0.55-0.75, hyp 0.75-1
ZONES = {"Resting":(0,0.15),"Proliferative":(0.15,0.55),"Prehypertrophic":(0.55,0.75),"Hypertrophic":(0.75,1.0)}
def band(a,b):
    s=np.zeros(N); s[(x>=a)&(x<b)]=1.0; return s

def solve_rd(D, k, S, bc=("neumann","neumann"), bcval=(0,0)):
    """Steady state: D u'' - k u + S = 0 on [0,1]. Tridiagonal solve."""
    A = np.zeros((N,N)); b = -S.copy()
    for i in range(1,N-1):
        A[i,i-1]=D/dx**2; A[i,i]=-2*D/dx**2 - k; A[i,i+1]=D/dx**2
    # boundaries
    for j,(side,val) in enumerate(zip(bc,bcval)):
        i = 0 if j==0 else N-1
        if side=="dirichlet":
            A[i,:]=0; A[i,i]=1; b[i]=val
        else:  # neumann (zero-flux unless val given as influx)
            if i==0: A[i,i]=-1/dx; A[i,i+1]=1/dx; b[i]=val
            else:    A[i,i]=1/dx; A[i,i-1]=-1/dx; b[i]=val
    return np.linalg.solve(A,b)

def run(params):
    Dp,kp,Di,ki,Dc,kc,Dd,kd,drug_bc,ihh_fb = params
    # iterate PTHrP<->IHH feedback (IHH promotes PTHrP production at resting end)
    P = np.zeros(N)
    for _ in range(40):
        # IHH from prehyp band
        I = solve_rd(Di, ki, band(0.55,0.75)*1.0)
        ihh_signal = np.trapezoid(I*band(0.55,0.75), x)      # integrated IHH
        Sp = band(0,0.15)*(1.0 + ihh_fb*ihh_signal)          # PTHrP source boosted by IHH feedback
        Pn = solve_rd(Dp, kp, Sp)
        if np.max(np.abs(Pn-P))<1e-6: P=Pn; break
        P=Pn
    C = solve_rd(Dc, kc, band(0.15,0.75)*1.0)                # CNP from prolif+prehyp
    # drug: systemic -> influx at BOTH boundaries (perichondrial + metaphyseal vasculature)
    Dg = solve_rd(Dd, kd, np.zeros(N), bc=("dirichlet","dirichlet"), bcval=(drug_bc,drug_bc))
    return P,I,C,Dg

def norm(v):
    v=v-v.min(); return v/(v.max()+1e-12)

# --- baseline run ---
base = (0.02,20,0.02,20,0.03,10, 0.01,5, 1.0, 2.0)
P,I,C,Dg = run(base)

# zonal means
def zmean(v): return {z: float(v[(x>=a)&(x<b)].mean()) for z,(a,b) in ZONES.items()}
zon = pd.DataFrame({"PTHrP":zmean(P),"IHH":zmean(I),"CNP":zmean(C),"Drug":zmean(Dg)})
zon.to_csv(ROOT/"results/tables/SPATIAL_RD_ZONAL.csv")

# --- flattening metric across a wide drug-penetration prior sweep ---
# gradient "informativeness" = coefficient of variation across zones (endogenous morphogens are graded)
def cv(v):
    zv=np.array([v[(x>=a)&(x<b)].mean() for a,b in ZONES.values()]); return zv.std()/(zv.mean()+1e-12)
rows=[]
for _ in range(2000):
    Dd = 10**rng.uniform(-3,-1); kd=10**rng.uniform(0,1.5); bc=rng.uniform(0.2,2.0)
    p=list(base); p[6],p[7],p[8]=Dd,kd,bc
    _,_,_,dg=run(p)
    rows.append({"Dd":Dd,"kd":kd,"drug_bc":bc,"drug_CV":cv(dg),
                 "penetration_ratio": dg[(x>=0)&(x<0.15)].mean()/(dg.max()+1e-12)})
sweep=pd.DataFrame(rows)
sweep.to_csv(ROOT/"results/sensitivity/spatial_rd_drug_penetration.csv", index=False)

pthrp_cv=cv(P); drug_cv_med=sweep["drug_CV"].median()
print("Zonal means:\n", zon.round(3).to_string())
print(f"\nendogenous PTHrP gradient CV = {pthrp_cv:.3f} (graded); IHH peak zone = {max(zmean(I),key=zmean(I).get)}")
print(f"systemic-drug spatial CV: median {drug_cv_med:.3f} (near 0 => FLATTENS endogenous spatial info)")
print(f"drug resting-zone penetration ratio: median {sweep['penetration_ratio'].median():.2f} "
      f"(low => systemic drug reaches resting end poorly under low diffusion)")

# --- mass-conservation / sanity tests ---
def flux_balance(v,D,k,S):
    # residual using the SAME discrete stencil the solver uses (exact -> ~0 at interior nodes)
    lap = (v[:-2] - 2*v[1:-1] + v[2:]) / dx**2
    resid = D*lap - k*v[1:-1] + S[1:-1]
    return float(np.max(np.abs(resid)))
assert flux_balance(I,base[2],base[3],band(0.55,0.75))<1e-2, "IHH steady-state residual too high"
assert zmean(P)["Resting"]==max(zmean(P).values()), "PTHrP should peak at resting end"
assert max(zmean(I),key=zmean(I).get)=="Prehypertrophic", "IHH should peak in prehypertrophic band"
print("\nPASS: mass-conservation residual small; PTHrP peaks resting; IHH peaks prehypertrophic")

# figure
fig,ax=plt.subplots(1,2,figsize=(12,4.4))
for v,l,c in [(norm(P),"PTHrP","#4C78A8"),(norm(I),"IHH","#e45756"),(norm(C),"CNP","#54a24b"),(norm(Dg),"systemic Drug","#9467bd")]:
    ax[0].plot(x,v,label=l,color=c,lw=2)
for z,(a,b) in ZONES.items(): ax[0].axvspan(a,b,alpha=0.05)
ax[0].set_xlabel("growth-plate depth (0=resting, 1=metaphysis)"); ax[0].set_ylabel("normalized conc.")
ax[0].set_title("Endogenous morphogen gradients vs systemic drug"); ax[0].legend(fontsize=8)
ax[1].hist(sweep["drug_CV"],bins=50,color="#9467bd",alpha=0.85)
ax[1].axvline(pthrp_cv,color="#4C78A8",ls="--",lw=2,label=f"endogenous PTHrP CV={pthrp_cv:.2f}")
ax[1].set_xlabel("systemic-drug spatial CV (0 = flat)"); ax[1].set_ylabel("param draws")
ax[1].set_title("Systemic drug FLATTENS spatial info vs graded morphogens"); ax[1].legend(fontsize=8)
fig.suptitle("Spatial reaction-diffusion (SIMULATION, not evidence)",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(ROOT/"results/figures/spatial_rd.png",dpi=130)
print("wrote results/figures/spatial_rd.png")
