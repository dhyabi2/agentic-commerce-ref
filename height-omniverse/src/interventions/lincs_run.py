#!/usr/bin/env python
"""Phase F2 (partial): build the EMPIRICAL growth-plate 'desired-state' up/down signature from the
GSE288529 zonal atlas (resting/proliferative = productive-plate state vs terminal hypertrophic).
This signature is the INPUT a LINCS/CMap connectivity query needs. The L1000 sig-matrix itself
(GSE92742/GSE70138, ~50 GB) is NOT ingested here (see docs/BLOCKED_DATA.md); this script writes the
signature + a reproducible connectivity stub. Connectivity = hypothesis generation only; NEVER efficacy.
Most L1000 cell lines are not growth-plate chondrocytes - cell-context mismatch must be flagged."""
import pandas as pd, numpy as np, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "tables"
spec = pd.read_parquet(OUT / "GENOME_WIDE_STATE_SPECIFICITY.parquet")

# "productive plate" desired state = high in resting+proliferative, low in terminal hypertrophic
spec["productive"] = spec["mean_Resting"] + spec["mean_Proliferative"]
spec["terminal"] = spec["mean_Hypertrophic"]
spec["log2fc_productive_vs_terminal"] = np.log2((spec["productive"]+1e-3)/(spec["terminal"]+1e-3))
expressed = spec[(spec["productive"]+spec["terminal"]) > 0.3].copy()

UP = expressed.sort_values("log2fc_productive_vs_terminal", ascending=False).head(150)   # resting/prolif markers
DOWN = expressed.sort_values("log2fc_productive_vs_terminal").head(150)                    # terminal-hypertrophic
sig = pd.concat([
    UP.assign(direction="up_in_desired_state")[["gene","log2fc_productive_vs_terminal","direction"]],
    DOWN.assign(direction="down_in_desired_state")[["gene","log2fc_productive_vs_terminal","direction"]],
])
sig["basis"] = "GSE288529 zonal atlas (single library, descriptive)"
sig.to_csv(OUT / "GROWTH_PLATE_DESIRED_STATE_SIGNATURE.csv", index=False)

# LINCS connectivity stub (documented, not executed - data too large)
stub = pd.DataFrame([{
    "step":"L1000 connectivity","status":"externally_blocked (sig-matrix ~50GB)",
    "resource":"GSE92742 / GSE70138 (Broad LINCS L1000)",
    "method":"weighted connectivity of this UP/DOWN signature vs L1000 perturbagen signatures",
    "caveat":"most L1000 lines are NOT growth-plate chondrocytes; cell-context distance must be scored; hypothesis-generation only",
    "reproducible_by":"download L1000 Level5 sig-matrix + cmapPy; score connectivity; consensus across cell/dose/time"}])
stub.to_csv(OUT / "LINCS_CONNECTIVITY_STATUS.csv", index=False)
print(f"desired-state signature: {len(UP)} up / {len(DOWN)} down genes")
print("top UP (resting/proliferative):", UP["gene"].head(12).tolist())
print("top DOWN (terminal hypertrophic):", DOWN["gene"].head(12).tolist())
print("L1000 connectivity: externally_blocked (sig-matrix ~50GB); signature written for future query")
