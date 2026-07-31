#!/usr/bin/env python
"""
Data-driven cluster -> state assignment, independent of Pth1r.
Each Leiden cluster is scored against every signature panel using the mean of
per-gene z-scores (z across cluster means). Cluster assigned to argmax panel.
Contaminant panels compete on equal footing so non-growth-plate clusters fall out.
"""
import numpy as np
import pandas as pd

means = pd.read_csv("tables/cluster_marker_means.csv", index_col=0)
means.index = means.index.astype(str)

panels = {
    "Resting":              ["Apoe", "Clu", "C1qtnf3", "Pthlh"],
    "Proliferative":        ["Acan", "Col2a1", "Sox9", "Ccnd1"],
    "Prehypertrophic":      ["Ihh", "Runx2", "Sp7"],
    "Hypertrophic":         ["Col10a1", "Mmp13", "Ibsp"],
    "Articular":            ["Prg4", "Gdf5", "Cytl1"],
    "Osteo/Perichondrial":  ["Col1a1", "Col3a1", "Postn", "Bglap", "Alpl"],
    "Immune":               ["Ptprc", "C1qa", "Cd68", "Lyz2"],
    "Endothelial":          ["Pecam1", "Cdh5", "Kdr"],
}

# z-score each gene across clusters
z = (means - means.mean(axis=0)) / (means.std(axis=0) + 1e-9)

scores = pd.DataFrame(index=means.index)
for name, genes in panels.items():
    g = [x for x in genes if x in z.columns]
    scores[name] = z[g].mean(axis=1)

assign = scores.idxmax(axis=1)
out = scores.round(2).copy()
out["ASSIGNED"] = assign
out["Pth1r_mean"] = means["Pth1r"].round(2)
out["n_top2_markers"] = [", ".join(scores.loc[c].sort_values(ascending=False).index[:2]) for c in scores.index]
pd.set_option("display.width", 240, "display.max_columns", 40)
print(out.to_string())
print("\nAssignment counts:")
print(assign.value_counts().to_string())
out.to_csv("tables/cluster_assignment_scores.csv")
