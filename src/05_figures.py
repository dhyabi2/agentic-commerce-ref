#!/usr/bin/env python
"""Figures for the HDAC4 growth plate analysis.

Palette taken from the dataviz reference instance and validated with
scripts/validate_palette.js:
  categorical (3 slots, all-pairs): #2a78d6 #eb6834 #1baf7a  -> ALL CHECKS PASS
  ordinal blue ramp (4 steps):      #86b6ef #3987e5 #1c5cab #0d366b -> ALL CHECKS PASS
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hdac4_config as C  # noqa: E402

warnings.filterwarnings("ignore")

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8a86"
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
ZONE_RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
ZONE_COL = dict(zip(C.ZONE_ORDER, ZONE_RAMP))

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": "#d8d7d2", "grid.color": "#ececE8", "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 9, "axes.titlesize": 10, "figure.dpi": 150, "lines.linewidth": 2,
})

CAVEAT = ("Inference from target-gene transcript abundance only. "
          "RNA-seq cannot measure HDAC4 nuclear/cytoplasmic localisation.")


def _cap(fig, text: str) -> None:
    fig.text(0.005, 0.004, text, fontsize=6.5, color=MUTED, ha="left", va="bottom")


def fig_zonal(pb: pd.DataFrame) -> None:
    tags = ["human_GSE288028", "mouse_GSE288028", "mouse_GSE288529"]
    methods = [("target_primary_scanpy", "scanpy score_genes"), ("target_primary_ulm", "decoupler ULM")]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.4), sharex=True)
    for r, (col, mname) in enumerate(methods):
        for c, tag in enumerate(tags):
            ax = axes[r, c]
            d = pb[(pb.dataset_tag == tag) & pb.zone.isin(C.ZONE_ORDER)]
            if d.empty or col not in d or d[col].isna().all():
                ax.text(0.5, 0.5, "no data", ha="center", va="center", color=MUTED,
                        transform=ax.transAxes)
                ax.set_xticks(range(4)); ax.set_xticklabels([z[:6] for z in C.ZONE_ORDER])
                continue
            for sample, g in d.groupby("sample_id", observed=True):
                g = g.set_index("zone").reindex(C.ZONE_ORDER)
                ok = g[col].notna() & ~g["below_min_cells_flag"].fillna(True).astype(bool)
                ax.plot(np.arange(4)[ok.to_numpy()], g[col][ok.to_numpy()],
                        color=CAT[0], alpha=0.35, lw=1.2, marker="o", ms=3.5,
                        markeredgecolor=SURFACE, markeredgewidth=0.6)
            m = d[~d.below_min_cells_flag].groupby("zone", observed=True)[col].mean().reindex(C.ZONE_ORDER)
            ax.plot(range(4), m, color=CAT[1], lw=2.6, marker="o", ms=7,
                    markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=5,
                    label="mean across samples")
            ax.grid(axis="y", alpha=0.7); ax.set_axisbelow(True)
            ax.set_xticks(range(4))
            ax.set_xticklabels(["rest", "prolif", "prehyp", "hyper"])
            if r == 0:
                ax.set_title(tag.replace("_", " "), color=INK)
            if c == 0:
                ax.set_ylabel(f"target expression\n({mname})")
            if r == 0 and c == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.suptitle("A. Zonal gradient of HDAC4 repression-target expression "
                 "(primary 5-gene set, disjoint from zone markers)", y=0.99, color=INK)
    fig.text(0.5, 0.945, "higher target expression = lower inferred HDAC4 repressive activity; "
             "thin lines = individual libraries", ha="center", fontsize=7.5, color=INK2)
    _cap(fig, CAVEAT)
    fig.tight_layout(rect=[0, 0.015, 1, 0.935])
    fig.savefig(C.FIGS / "fig1_zonal_gradient.png", bbox_inches="tight")
    plt.close(fig)


def fig_score_vs_age(pb: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))

    ax = axes[0]
    d = pb[(pb.dataset_tag == "human_GSE288028") & (pb.zone == "resting")
           & (pb.arm == "uncultured")]
    if not d.empty:
        donors = sorted(d.donor.unique())
        for i, dn in enumerate(donors):
            g = d[d.donor == dn]
            flag = bool(g.below_min_cells_flag.iloc[0])
            ax.scatter([i], g["target_primary_scanpy"], s=70, color=CAT[0] if not flag else MUTED,
                       edgecolor=SURFACE, linewidth=1.2, zorder=4)
            ax.annotate(f"n={int(g.n_cells.iloc[0])}", (i, g["target_primary_scanpy"].iloc[0]),
                        textcoords="offset points", xytext=(0, 11), ha="center",
                        fontsize=7, color=INK2)
        ax.set_xticks(range(len(donors)))
        ax.set_xticklabels(donors, rotation=20, ha="right", fontsize=7.5)
    ax.set_ylabel("resting-zone target expression\n(primary set, scanpy)")
    ax.set_xlabel("donor  —  NOT an age axis")
    ax.set_title("B. Human: no age variable exists", color=INK)
    ax.grid(axis="y", alpha=0.7); ax.set_axisbelow(True)
    ax.margins(y=0.28)
    ax.text(0.5, 0.03, "GEO records no per-donor age or Tanner stage.\n"
                       "Donor order is NOT a proxy for age. No age model fitted.\n"
                       "grey = zone flagged <100 cells",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7.5, color="#b4451f",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fdf1ec", ec="#eb6834", lw=1))

    ax = axes[1]
    p = C.TABLES / "E_rat_sample_scores.tsv"
    if p.exists():
        r = pd.read_csv(p, sep="\t").dropna(subset=["age_weeks", "target_primary_zmean"])
        # the 3-12 wk animals give one combined proliferative+early-hypertrophic
        # fraction; the 1 wk animals were dissected into separate zones
        zones = list(dict.fromkeys(r.sort_values("age_weeks").zone))
        # distinct step per series, in age/differentiation order (never two the same)
        colmap = {z: ZONE_RAMP[i % len(ZONE_RAMP)] for i, z in enumerate(zones)}
        for i, z in enumerate(zones):
            g = r[r.zone == z].sort_values("age_weeks")
            jit = (i - (len(zones) - 1) / 2) * 0.18
            ax.scatter(g.age_weeks + jit, g.target_primary_zmean, s=34,
                       color=colmap[z], edgecolor=SURFACE, linewidth=0.8,
                       label=z.replace(" zone", "").replace(
                           "proliferative and early hypertrophic", "prolif + early hyper"),
                       zorder=4)
            m = g.groupby("age_weeks")["target_primary_zmean"].mean()
            if len(m) > 1:
                ax.plot(m.index + jit, m.values, color=colmap[z], lw=2, alpha=0.85)
        # fitted age slope on the castrated series
        e = pd.read_csv(C.TABLES / "E_rat_age_effect.tsv", sep="\t")
        e = e[(e.score == "target_primary_zmean") & (e.subset == "castrated_only")]
        if len(e):
            row = e.iloc[0]
            xs = np.array([3, 12], dtype=float)
            sub = r[r.castrated == "yes"]
            b = sub.target_primary_zmean.mean() - row.slope_per_week * sub.age_weeks.mean()
            ax.plot(xs, row.slope_per_week * xs + b, color=CAT[1], lw=2.4, ls="--",
                    zorder=5, label=f"fit 3–12 wk: {row.slope_per_week:+.3f}/wk (p={row.p_value:.4f})")
        ax.legend(frameon=False, fontsize=7, title="zone", title_fontsize=7.5, loc="upper left")
        ax.set_xlabel("age (weeks)   —   1 wk intact, 3–12 wk castrated at day 18")
        ax.set_ylabel("target expression (z-mean, primary set)")
        ax.set_title("Rat GSE16981 (supplementary): the only real age series", color=INK)
        ax.grid(axis="y", alpha=0.7); ax.set_axisbelow(True)
    else:
        ax.text(0.5, 0.5, "rat analysis not run", ha="center", va="center", color=MUTED)
    _cap(fig, CAVEAT + "  Rat is not a fusing species; do not extrapolate to human fusion.")
    fig.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(C.FIGS / "fig2_score_vs_age.png", bbox_inches="tight")
    plt.close(fig)


def fig_boundary(cross: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.6))
    d = cross[(cross.get("status") == "ok") & (cross.gene_set == "primary")
              & (cross.method == "scanpy")] if "status" in cross else pd.DataFrame()
    if d.empty:
        ax.text(0.5, 0.5, "no sample yielded a boundary crossing\n(see C_boundary_crossing.tsv)",
                ha="center", va="center", color=MUTED)
    else:
        d = d.sort_values(["dataset_tag", "donor", "arm"])
        arms = list(dict.fromkeys(d.arm))
        cmap = {a: CAT[i % 3] for i, a in enumerate(arms)}
        y = np.arange(len(d))
        ax.hlines(y, 0, d.crossing_pseudotime_quantile, color="#e2e1dc", lw=2)
        for a in arms:
            m = (d.arm == a).to_numpy()
            ax.scatter(d.crossing_pseudotime_quantile[m], y[m], s=64, color=cmap[a],
                       edgecolor=SURFACE, linewidth=1.2, label=a, zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.donor} · {r.arm}" for r in d.itertuples()], fontsize=7.5)
        ax.set_xlim(0, 1)
        ax.set_xlabel("crossing position along differentiation pseudotime "
                      "(0 = resting end, 1 = hypertrophic end)")
        ax.legend(frameon=False, fontsize=8, title="arm", title_fontsize=8)
        ax.grid(axis="x", alpha=0.7); ax.set_axisbelow(True)
    ax.set_title("C. Boundary crossing position per library\n"
                 "(cannot be tested against age — no age variable exists)", color=INK)
    _cap(fig, CAVEAT)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(C.FIGS / "fig3_boundary_crossing.png", bbox_inches="tight")
    plt.close(fig)


def fig_pathway(pb: pd.DataFrame) -> None:
    tags = ["human_GSE288028", "mouse_GSE288028", "mouse_GSE288529"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.2))
    for ax, tag in zip(axes, tags):
        d = pb[(pb.dataset_tag == tag) & pb.zone.isin(C.ZONE_ORDER) & ~pb.below_min_cells_flag]
        cols = [c for c in d.columns if c.startswith("expr_")]
        if d.empty or not cols:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", color=MUTED); ax.axis("off")
            continue
        m = d.groupby("zone", observed=True)[cols].mean().reindex(C.ZONE_ORDER).dropna(how="all")
        z = (m - m.mean()) / m.std(ddof=0).replace(0, np.nan)
        z = z.dropna(axis=1, how="all")
        names = [c.replace("expr_", "") for c in z.columns]
        im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=90, fontsize=7)
        ax.set_yticks(range(len(z.index)))
        ax.set_yticklabels(list(z.index), fontsize=8)
        ax.set_title(tag.replace("_", " "), fontsize=9, color=INK)
        for s in ax.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02)
    cb.set_label("z-score of mean expression across zones", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("Pathway-component expression by zone (reported directly, never folded "
                 "into the activity score)", color=INK, y=1.0)
    _cap(fig, CAVEAT)
    fig.savefig(C.FIGS / "fig4_pathway_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def fig_controls(pb: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))

    ax = axes[0]
    d = pb[(pb.dataset_tag == "human_GSE288028") & pb.zone.isin(C.ZONE_ORDER)
           & ~pb.below_min_cells_flag]
    for i, (col, lab) in enumerate([("target_primary_scanpy", "HDAC4 target set (primary)"),
                                    ("target_housekeeping_scanpy", "housekeeping control set")]):
        if col not in d:
            continue
        m = d.groupby("zone", observed=True)[col].agg(["mean", "sem"]).reindex(C.ZONE_ORDER)
        ax.errorbar(range(4), m["mean"], yerr=m["sem"], color=CAT[i], lw=2, marker="o",
                    ms=7, capsize=3, markeredgecolor=SURFACE, markeredgewidth=1.2, label=lab)
    ax.set_xticks(range(4)); ax.set_xticklabels(["rest", "prolif", "prehyp", "hyper"])
    ax.set_ylabel("score"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("Negative control: housekeeping set must be flat", color=INK)
    ax.grid(axis="y", alpha=0.7); ax.set_axisbelow(True)

    ax = axes[1]
    p = C.TABLES / "ctrl_permutation_zone_shuffle.tsv"
    if p.exists():
        c = pd.read_csv(p, sep="\t")
        row = c[(c.gene_set == "primary") & (c.method == "scanpy")]
        if len(row):
            r = row.iloc[0]
            ax.axvspan(r["null_p2.5"], r["null_p97.5"], color="#e8eef8", label="null 95% band")
            ax.axvline(r["null_mean"], color=MUTED, lw=2, label="null mean (labels shuffled)")
            ax.axvline(r["observed_mean_rho"], color=CAT[1], lw=2.6, label="observed")
            ax.set_xlabel("mean Spearman rho (zone rank vs score)")
            ax.set_yticks([])
            ax.set_title(f"Zone-label shuffle, {int(r['n_permutations'])} permutations\n"
                         f"empirical p = {r['empirical_p_two_sided']:.4f}", color=INK)
            ax.legend(frameon=False, fontsize=7.5)
    else:
        ax.text(0.5, 0.5, "permutation control not run", ha="center", va="center", color=MUTED)
    _cap(fig, CAVEAT)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(C.FIGS / "fig5_controls.png", bbox_inches="tight")
    plt.close(fig)


def fig_umap() -> None:
    import scanpy as sc
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
    for ax, tag in zip(axes, ["human_GSE288028", "mouse_GSE288028", "mouse_GSE288529"]):
        p = C.PROC / f"{tag}.scored.h5ad"
        if not p.exists():
            ax.axis("off"); continue
        ad = sc.read_h5ad(p)
        if "X_umap" not in ad.obsm:
            ax.axis("off"); continue
        u = ad.obsm["X_umap"]
        order = [*C.ZONE_ORDER, "non_chondrocyte"]
        for z in order:
            m = (ad.obs.zone == z).to_numpy()
            if not m.any():
                continue
            ax.scatter(u[m, 0], u[m, 1], s=1.2, linewidths=0,
                       color=ZONE_COL.get(z, MUTED), label=f"{z} ({int(m.sum())})", alpha=0.7)
        ax.set_title(tag.replace("_", " "), fontsize=9, color=INK)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.legend(frameon=False, fontsize=6.5, markerscale=6, loc="best")
        del ad
    fig.suptitle("Zone assignment (from zone markers only, independent of the scoring genes)",
                 color=INK)
    _cap(fig, CAVEAT)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(C.FIGS / "fig6_umap_zones.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    pb = pd.read_csv(C.TABLES / "pseudobulk_scores.tsv", sep="\t")
    fig_zonal(pb); print("fig1 ok", flush=True)
    fig_score_vs_age(pb); print("fig2 ok", flush=True)
    cp = C.TABLES / "C_boundary_crossing.tsv"
    fig_boundary(pd.read_csv(cp, sep="\t") if cp.exists() else pd.DataFrame())
    print("fig3 ok", flush=True)
    fig_pathway(pb); print("fig4 ok", flush=True)
    fig_controls(pb); print("fig5 ok", flush=True)
    fig_umap(); print("fig6 ok", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
