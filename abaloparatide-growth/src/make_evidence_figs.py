#!/usr/bin/env python
"""Evidence-map + directness-hierarchy figures from the ledger."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

el = pd.read_csv("results/tables/EVIDENCE_LEDGER.csv")
tier_order = ["T1","T2","T3","T4","T5","T6","T7","T8"]
tier_label = {
 "T1":"T1 abaloparatide | conv. physis | FINAL length",
 "T2":"T2 abaloparatide | conv. physis | short-term",
 "T3":"T3 abaloparatide | mandibular condyle",
 "T4":"T4 abaloparatide | fracture/repair",
 "T5":"T5 other PTH1R agonist | normal juvenile physis",
 "T6":"T6 PTH1R agonist | disease/replacement (peds)",
 "T7":"T7 human genetics / transgenic (chronic)",
 "T8":"T8 molecular / cellular mechanism",
}
counts = {t: int((el["tier"]==t).sum()) for t in tier_order}
# directness score (closeness to primary estimand): T1 highest ... T8 lowest
directness = {"T1":1.0,"T2":0.8,"T3":0.35,"T4":0.3,"T5":0.55,"T6":0.5,"T7":0.35,"T8":0.25}

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(15,6.2), gridspec_kw={"width_ratios":[1.25,1]})

# --- left: directness ladder ---
y = np.arange(len(tier_order))[::-1]
for ti,t in zip(y, tier_order):
    n = counts[t]; d = directness[t]
    color = "#d62728" if (t in ("T1","T2") and n==0) else plt.cm.viridis(d)
    ax1.barh(ti, d, color=color, edgecolor="k", alpha=0.9)
    tag = f"  n={n}" + ("  <-- NO DIRECT EVIDENCE" if (t in ("T1","T2") and n==0) else f"  n={n}") if False else f"  n={n}"
    lbl = tier_label[t] + ("   *** EMPTY ***" if (t in ("T1","T2") and n==0) else "")
    ax1.text(0.01, ti, lbl, va="center", ha="left", fontsize=8.5, color="white" if d>0.5 else "black")
    ax1.text(d+0.01, ti, f"n={n}", va="center", ha="left", fontsize=9, fontweight="bold")
ax1.set_yticks([]); ax1.set_xlim(0,1.25); ax1.set_xlabel("directness to primary estimand (final conventional long-bone length)")
ax1.set_title("Evidence directness hierarchy\nThe two tiers that would answer the question directly are EMPTY", fontsize=11)

# --- right: comparator calibration (velocity -> final height translation) ---
comp = ["GH (ISS)","Vosoritide\n(achondro.)","Letrozole\n(AI, boys)","Abaloparatide\n(this analysis)"]
velocity = [1,1,1,0.5]     # qualitative: did it raise velocity/predicted?
finalht  = [1,np.nan,0,np.nan]  # did it raise ATTAINED final height? nan=unknown
x = np.arange(len(comp))
ax2.bar(x-0.18, velocity, 0.36, label="raised velocity / predicted", color="#4C78A8")
fh_plot = [v if not np.isnan(v) else 0 for v in finalht]
hatch = ["","//","","//"]
bars = ax2.bar(x+0.18, fh_plot, 0.36, label="raised ATTAINED final height", color="#54a24b")
for i,v in enumerate(finalht):
    if np.isnan(v): ax2.text(x[i]+0.18, 0.05, "?", ha="center", fontsize=13, fontweight="bold", color="k")
    elif v==0: ax2.text(x[i]+0.18, 0.03, "NO", ha="center", fontsize=9, fontweight="bold", color="#d62728")
ax2.set_xticks(x); ax2.set_xticklabels(comp, fontsize=8.5); ax2.set_ylim(0,1.2)
ax2.set_yticks([0,1]); ax2.set_yticklabels(["no","yes"])
ax2.legend(fontsize=8, loc="upper right")
ax2.set_title("Calibration: velocity/predicted gains do NOT reliably\nbecome attained final height (letrozole = cautionary case)", fontsize=10.5)

fig.suptitle("Abaloparatide & final long-bone length — evidence landscape", fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig("results/figures/evidence_map.png", dpi=130)
print("wrote results/figures/evidence_map.png; tier counts:", counts)
