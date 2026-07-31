#!/usr/bin/env python
"""Primary analyses A-D and the pre-registered negative controls.

A  Zonal gradient (sanity check)      -> results/tables/A_zonal_gradient*.tsv
B  Age effect                         -> results/tables/B_age_effect*.tsv
C  Boundary migration                 -> results/tables/C_boundary_crossing.tsv
D  GH arm                             -> results/tables/D_gh_contrast*.tsv
6  Negative controls and guards       -> results/tables/ctrl_*.tsv

Sign convention: all "target_*" scores measure TARGET EXPRESSION. Inferred HDAC4
repressive activity is the negative. The zonal prediction is that target expression
increases resting -> hypertrophic.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hdac4_config as C  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.n_jobs = 4
RNG = np.random.default_rng(C.SEED)

METHODS = ["scanpy", "ulm"]
SETNAMES = ["primary", "secondary", "housekeeping"]
ZONE_RANK = {z: i for i, z in enumerate(C.ZONE_ORDER)}

RESULTS: dict = {}


# ======================================================================================
# A. Zonal gradient
# ======================================================================================
def zonal_gradient(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Spearman rho between zone rank and pseudobulk mean score, per sample.

    Run twice: with the pre-registered <100-cell zone filter applied, and without it.
    The human libraries contain almost no hypertrophic chondrocytes, so the filtered
    pass drops most of the zonal dynamic range; both passes are reported so the reader
    can see what the filter costs.
    """
    rows = []
    for filt in (True, False):
        d = pb[pb.zone.isin(C.ZONE_ORDER)].copy()
        if filt:
            d = d[~d.below_min_cells_flag]
        d["zone_rank"] = d.zone.map(ZONE_RANK)

        for (tag, sample, donor, arm), g in d.groupby(
            ["dataset_tag", "sample_id", "donor", "arm"], observed=True
        ):
            if g.zone.nunique() < 3:
                continue
            for sn in SETNAMES:
                for m in METHODS:
                    col = f"target_{sn}_{m}"
                    if col not in g or g[col].isna().all():
                        continue
                    rho, p = stats.spearmanr(g.zone_rank, g[col])
                    lo = g.loc[g.zone == "resting", col]
                    hi = g.loc[g.zone == "hypertrophic", col]
                    rows.append(
                        {
                            "dataset_tag": tag, "sample_id": sample, "donor": donor,
                            "arm": arm, "gene_set": sn, "method": m,
                            "min_cell_filter_applied": filt,
                            "n_zones": g.zone.nunique(),
                            "spearman_rho": rho, "spearman_p": p,
                            "resting_mean": lo.mean() if len(lo) else np.nan,
                            "hypertrophic_mean": hi.mean() if len(hi) else np.nan,
                            "hyper_minus_resting": (hi.mean() - lo.mean())
                            if len(lo) and len(hi) else np.nan,
                        }
                    )
    per_sample = pd.DataFrame(rows)

    # aggregate to the biological replicate (donor), then test across donors
    agg = []
    for (tag, sn, m, filt), g in per_sample.groupby(
        ["dataset_tag", "gene_set", "method", "min_cell_filter_applied"]
    ):
        per_donor = g.groupby("donor")["spearman_rho"].mean()
        n = len(per_donor)
        mean_rho = per_donor.mean()
        if n >= 2 and per_donor.std(ddof=1) > 0:
            t = stats.ttest_1samp(per_donor, 0.0)
            se = per_donor.std(ddof=1) / np.sqrt(n)
            tcrit = stats.t.ppf(0.975, n - 1)
            ci = (mean_rho - tcrit * se, mean_rho + tcrit * se)
            tstat, pval = t.statistic, t.pvalue
        else:
            ci, tstat, pval = (np.nan, np.nan), np.nan, np.nan
        n_pos = int((per_donor > 0).sum())
        sign_p = stats.binomtest(n_pos, n, 0.5).pvalue if n else np.nan
        agg.append(
            {
                "dataset_tag": tag, "gene_set": sn, "method": m,
                "min_cell_filter_applied": filt,
                "n_donors": n, "mean_rho_across_donors": mean_rho,
                "ci95_lo": ci[0], "ci95_hi": ci[1],
                "t_stat": tstat, "t_p": pval,
                "n_donors_positive": n_pos, "sign_test_p": sign_p,
                "per_donor_rho": json.dumps({k: round(float(v), 4) for k, v in per_donor.items()}),
            }
        )
    return per_sample, pd.DataFrame(agg)


# ======================================================================================
# B. Age effect -- adequacy gate
# ======================================================================================
def age_effect(pb: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    ages = {d: a for d, a in C.DONOR_AGE_YEARS.items() if a is not None}
    gate = {
        "gate_a_per_donor_age_recorded": bool(ages),
        "n_donors_with_age": len(ages),
        "series_level_age_statement": C.AGE_RANGE_SERIES_LEVEL,
        "gate_b_age_span_adequate": False,
        "verdict": None,
    }
    if not ages:
        gate["verdict"] = "UNDEFINED - regressor does not exist"
        gate["note"] = (
            "No per-donor age or Tanner stage is recorded in GEO for GSE288028 "
            "(verified in SOFT and MINiML). Per the pre-registration, no age model is "
            "fitted and no age surrogate (donor index, accession order, batch order) "
            "is substituted."
        )

    # What CAN be quantified: the donor-to-donor spread the RZ score actually shows.
    # Any real age effect would have to exceed this to be detectable at n=4.
    rz = pb[
        (pb.zone == "resting")
        & (pb.arm == "uncultured")
        & (pb.dataset_tag == "human_GSE288028")
        & (~pb.below_min_cells_flag)
    ]
    rows = []
    for sn in SETNAMES:
        for m in METHODS:
            col = f"target_{sn}_{m}"
            if col not in rz or rz[col].isna().all():
                continue
            v = rz.groupby("donor")[col].mean().dropna()
            n = len(v)
            sd = v.std(ddof=1) if n > 1 else np.nan
            rows.append(
                {
                    "gene_set": sn, "method": m, "n_donors": n,
                    "rz_mean": v.mean(), "rz_sd_between_donors": sd,
                    "rz_min": v.min(), "rz_max": v.max(),
                    "rz_range": v.max() - v.min(),
                    # half-width of a 95% CI on the donor mean at this n
                    "ci95_halfwidth_on_mean": (
                        stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n) if n > 1 else np.nan
                    ),
                    # smallest slope-per-year detectable if the (unavailable) ages had
                    # spanned the series-level 11-14y range, 80% power, alpha=.05
                    "min_detectable_diff_80pct_power": (
                        (stats.t.ppf(0.975, n - 2) + stats.t.ppf(0.80, n - 2))
                        * sd * np.sqrt(2 / n)
                        if n > 2 else np.nan
                    ),
                    "per_donor": json.dumps({k: round(float(x), 4) for k, x in v.items()}),
                }
            )
    return pd.DataFrame(rows), gate


def technical_bootstrap_rz(scored: dict[str, sc.AnnData], n_boot: int = 1000) -> pd.DataFrame:
    """Per-donor resting-zone mean with a bootstrap-over-cells interval.

    THIS IS A TECHNICAL INTERVAL. It bootstraps cells within one library and says
    nothing about donor-to-donor variability. Labelled as such in every column name.
    """
    ad = scored["human_GSE288028"]
    obs = ad.obs
    rows = []
    sel = obs[(obs.zone == "resting") & (obs.arm == "uncultured")]
    for donor, g in sel.groupby("donor", observed=True):
        for sn in SETNAMES:
            for m in METHODS:
                col = f"target_{sn}_{m}"
                if col not in g or g[col].isna().all():
                    continue
                x = g[col].to_numpy(dtype=float)
                x = x[~np.isnan(x)]
                if len(x) < 10:
                    continue
                idx = RNG.integers(0, len(x), size=(n_boot, len(x)))
                bs = x[idx].mean(1)
                rows.append(
                    {
                        "donor": donor, "gene_set": sn, "method": m,
                        "n_cells": len(x), "mean": x.mean(),
                        "TECHNICAL_ci95_lo": np.percentile(bs, 2.5),
                        "TECHNICAL_ci95_hi": np.percentile(bs, 97.5),
                        "interval_type": "TECHNICAL (bootstrap of cells within one library)",
                    }
                )
    return pd.DataFrame(rows)


# ======================================================================================
# C. Boundary migration
# ======================================================================================
def boundary_crossing(scored: dict[str, sc.AnnData]) -> pd.DataFrame:
    rows = []
    for tag, ad in scored.items():
        for sample, idx in ad.obs.groupby("sample_id", observed=True).groups.items():
            sub = ad[list(idx)].copy()
            if sub.n_obs < 500:
                rows.append({"dataset_tag": tag, "sample_id": sample,
                             "status": f"skipped: only {sub.n_obs} cells"})
                continue
            zones_present = set(sub.obs.zone) & set(C.ZONE_ORDER)
            if not {"resting", "hypertrophic"} <= zones_present:
                rows.append({"dataset_tag": tag, "sample_id": sample,
                             "status": "skipped: resting and/or hypertrophic zone absent"})
                continue
            try:
                sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca", random_state=C.SEED)
                sc.tl.diffmap(sub)
                root = int(np.argmax(sub.obs["zscore_resting"].to_numpy()))
                sub.uns["iroot"] = root
                sc.tl.dpt(sub)
                pt = sub.obs["dpt_pseudotime"].to_numpy(dtype=float)
            except Exception as exc:  # noqa: BLE001
                rows.append({"dataset_tag": tag, "sample_id": sample,
                             "status": f"skipped: pseudotime failed ({exc})"})
                continue
            ok = np.isfinite(pt)
            if ok.sum() < 200:
                rows.append({"dataset_tag": tag, "sample_id": sample,
                             "status": "skipped: pseudotime not finite"})
                continue
            # orient so that pseudotime increases resting -> hypertrophic
            zr = sub.obs.zone.map(ZONE_RANK).to_numpy(dtype=float)
            fin = ok & np.isfinite(zr)
            if fin.sum() > 10 and stats.spearmanr(pt[fin], zr[fin]).statistic < 0:
                pt = np.nanmax(pt[ok]) - pt
            rank = pd.Series(pt).rank(pct=True).to_numpy()

            for sn in ("primary", "secondary"):
                for m in METHODS:
                    col = f"target_{sn}_{m}"
                    if col not in sub.obs or sub.obs[col].isna().all():
                        continue
                    y = sub.obs[col].to_numpy(dtype=float)
                    good = ok & np.isfinite(y)
                    lo = np.nanmean(y[(sub.obs.zone == "resting").to_numpy() & good])
                    hi = np.nanmean(y[(sub.obs.zone == "hypertrophic").to_numpy() & good])
                    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
                        rows.append({
                            "dataset_tag": tag, "sample_id": sample, "gene_set": sn,
                            "method": m, "status": "no crossing: hypertrophic <= resting",
                            "resting_mean": lo, "hypertrophic_mean": hi})
                        continue
                    mid = 0.5 * (lo + hi)
                    bins = np.linspace(0, 1, 51)
                    which = np.clip(np.digitize(rank, bins) - 1, 0, 49)
                    prof = pd.Series(y[good]).groupby(which[good]).mean().reindex(range(50))
                    prof = prof.interpolate(limit_direction="both")
                    above = np.where(prof.to_numpy() >= mid)[0]
                    cross = (above[0] + 0.5) / 50.0 if len(above) else np.nan
                    rows.append(
                        {
                            "dataset_tag": tag, "sample_id": sample,
                            "donor": sub.obs.donor.iloc[0], "arm": sub.obs.arm.iloc[0],
                            "gene_set": sn, "method": m, "status": "ok",
                            "n_cells": int(good.sum()),
                            "resting_mean": lo, "hypertrophic_mean": hi, "midpoint": mid,
                            "crossing_pseudotime_quantile": cross,
                        }
                    )
    return pd.DataFrame(rows)


# ======================================================================================
# D. GH arm
# ======================================================================================
def identify_gp_subclusters(ad: sc.AnnData) -> pd.DataFrame:
    """Subcluster resting/progenitor cells and find which subcluster culture removes."""
    rz = ad[(ad.obs.zone == "resting")].copy()
    if rz.n_obs < 200:
        return pd.DataFrame()
    sc.pp.neighbors(rz, n_neighbors=15, use_rep="X_pca", random_state=C.SEED)
    sc.tl.leiden(rz, resolution=0.3, key_added="rz_sub", random_state=C.SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    rz.obs["cultured"] = rz.obs.arm.isin(["vehicle", "gh"])
    tab = (
        rz.obs.groupby(["rz_sub", "cultured"], observed=True).size()
        .unstack(fill_value=0)
    )
    frac = tab / tab.sum(0)
    out = pd.DataFrame(
        {
            "rz_sub": tab.index,
            "n_uncultured": tab.get(False, pd.Series(0, index=tab.index)).values,
            "n_cultured": tab.get(True, pd.Series(0, index=tab.index)).values,
            "frac_of_uncultured": frac.get(False, pd.Series(0, index=tab.index)).values,
            "frac_of_cultured": frac.get(True, pd.Series(0, index=tab.index)).values,
        }
    )
    out["culture_depletion_log2fc"] = np.log2(
        (out.frac_of_cultured + 1e-6) / (out.frac_of_uncultured + 1e-6)
    )
    cyc = rz.obs.groupby("rz_sub", observed=True)[["S_score", "G2M_score"]].mean()
    out = out.merge(cyc.reset_index(), on="rz_sub", how="left")
    out["cycling_score"] = out.S_score + out.G2M_score
    # GP1 = quiescent cluster most depleted by culture
    out = out.sort_values("culture_depletion_log2fc")
    out["label"] = ["GP2"] * len(out)
    if len(out) >= 2:
        cand = out[out.culture_depletion_log2fc < 0]
        if len(cand):
            gp1 = cand.sort_values("cycling_score").index[0]
            out.loc[gp1, "label"] = "GP1_quiescent_culture_depleted"
    # map back
    ad.obs["rz_sub"] = pd.Series(rz.obs["rz_sub"], index=rz.obs_names).reindex(ad.obs_names)
    lab = dict(zip(out.rz_sub, out.label))
    ad.obs["rz_gp_label"] = ad.obs["rz_sub"].map(lab)
    return out


def gh_contrast(ad: sc.AnnData, restrict_gp2: bool) -> pd.DataFrame:
    obs = ad.obs
    sel = obs[(obs.zone == "resting") & (obs.arm.isin(["vehicle", "gh"]))]
    if restrict_gp2:
        sel = sel[sel.rz_gp_label == "GP2"]
    rows = []
    for sn in SETNAMES:
        for m in METHODS:
            col = f"target_{sn}_{m}"
            if col not in sel or sel[col].isna().all():
                continue
            pbd = sel.groupby(["donor", "arm"], observed=True)[col].mean().unstack()
            pbd = pbd.dropna()
            if len(pbd) < 2 or "gh" not in pbd or "vehicle" not in pbd:
                continue
            diff = (pbd["gh"] - pbd["vehicle"]).astype(float)
            n = len(diff)
            sd = diff.std(ddof=1)
            se = sd / np.sqrt(n) if n > 1 else np.nan
            tcrit = stats.t.ppf(0.975, n - 1) if n > 1 else np.nan
            tt = stats.ttest_rel(pbd["gh"], pbd["vehicle"]) if n > 1 else None
            rows.append(
                {
                    "gene_set": sn, "method": m,
                    "restricted_to_GP2": restrict_gp2,
                    "n_donors_paired": n,
                    "mean_gh_minus_vehicle_target": diff.mean(),
                    "mean_gh_minus_vehicle_hdac4_activity": -diff.mean(),
                    "sd": sd,
                    "ci95_lo": diff.mean() - tcrit * se if n > 1 else np.nan,
                    "ci95_hi": diff.mean() + tcrit * se if n > 1 else np.nan,
                    "paired_t": tt.statistic if tt else np.nan,
                    "p_value": tt.pvalue if tt else np.nan,
                    "per_donor_diff": json.dumps(
                        {k: round(float(v), 4) for k, v in diff.items()}
                    ),
                }
            )
    return pd.DataFrame(rows)


# ======================================================================================
# 6. Negative controls
# ======================================================================================
def permutation_control(scored: dict[str, sc.AnnData], n_perm: int) -> pd.DataFrame:
    """Shuffle zone labels within each sample; the zonal gradient must vanish."""
    ad = scored["human_GSE288028"]
    obs = ad.obs
    keep = obs.zone.isin(C.ZONE_ORDER)
    rows = []
    for sn in ("primary", "secondary", "housekeeping"):
        for m in METHODS:
            col = f"target_{sn}_{m}"
            if col not in obs or obs[col].isna().all():
                continue
            observed_rhos, null = [], []
            per_sample_arrays = []
            for _, g in obs[keep].groupby("sample_id", observed=True):
                y = g[col].to_numpy(dtype=float)
                z = g.zone.map(ZONE_RANK).to_numpy()
                ok = np.isfinite(y)
                if ok.sum() < 200 or len(np.unique(z[ok])) < 3:
                    continue
                per_sample_arrays.append((y[ok], z[ok]))
                means = pd.Series(y[ok]).groupby(z[ok]).mean()
                observed_rhos.append(stats.spearmanr(means.index, means.values).statistic)
            if not per_sample_arrays:
                continue
            obs_mean = float(np.mean(observed_rhos))
            for _ in range(n_perm):
                rs = []
                for y, z in per_sample_arrays:
                    zp = RNG.permutation(z)
                    means = pd.Series(y).groupby(zp).mean()
                    rs.append(stats.spearmanr(means.index, means.values).statistic)
                null.append(np.mean(rs))
            null = np.asarray(null, dtype=float)
            rows.append(
                {
                    "gene_set": sn, "method": m, "n_samples": len(per_sample_arrays),
                    "observed_mean_rho": obs_mean,
                    "null_mean": float(np.nanmean(null)),
                    "null_sd": float(np.nanstd(null)),
                    "null_p2.5": float(np.nanpercentile(null, 2.5)),
                    "null_p97.5": float(np.nanpercentile(null, 97.5)),
                    "empirical_p_two_sided": float(
                        (np.sum(np.abs(null) >= abs(obs_mean)) + 1) / (len(null) + 1)
                    ),
                    "n_permutations": n_perm,
                }
            )
    return pd.DataFrame(rows)


def technical_covariate_control(pb: pd.DataFrame) -> pd.DataFrame:
    d = pb[(pb.dataset_tag == "human_GSE288028") & pb.zone.isin(C.ZONE_ORDER)
           & ~pb.below_min_cells_flag]
    rows = []
    for sn in SETNAMES:
        for m in METHODS:
            col = f"target_{sn}_{m}"
            if col not in d or d[col].isna().all():
                continue
            for cov in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
                r = stats.spearmanr(d[col], d[cov], nan_policy="omit")
                rows.append({"gene_set": sn, "method": m, "covariate": cov,
                             "spearman_rho": r.statistic, "p_value": r.pvalue,
                             "n_pseudobulk_units": len(d)})
    return pd.DataFrame(rows)


def leave_one_donor_out(pb: pd.DataFrame) -> pd.DataFrame:
    """Refit the zonal gradient leaving out each donor in turn.

    Run with and without the <100-cell zone filter. With the filter applied only one
    human donor retains three usable zones, so that pass carries no leave-one-out
    information at all -- which is itself worth showing.
    """
    rows = []
    for filt in (True, False):
        d = pb[(pb.dataset_tag == "human_GSE288028") & pb.zone.isin(C.ZONE_ORDER)].copy()
        if filt:
            d = d[~d.below_min_cells_flag]
        d["zone_rank"] = d.zone.map(ZONE_RANK)
        donors = sorted(pb[pb.dataset_tag == "human_GSE288028"].donor.unique())
        for sn in SETNAMES:
            for m in METHODS:
                col = f"target_{sn}_{m}"
                if col not in d or d[col].isna().all():
                    continue
                for drop in [None, *donors]:
                    sub = d if drop is None else d[d.donor != drop]
                    rr, dn = [], set()
                    for (sid, dnr), g in sub.groupby(["sample_id", "donor"], observed=True):
                        if g.zone.nunique() < 3:
                            continue
                        rr.append(stats.spearmanr(g.zone_rank, g[col]).statistic)
                        dn.add(dnr)
                    rows.append({"gene_set": sn, "method": m,
                                 "min_cell_filter_applied": filt,
                                 "donor_left_out": drop or "none (full)",
                                 "mean_rho": float(np.nanmean(rr)) if rr else np.nan,
                                 "n_samples": len(rr), "n_donors_contributing": len(dn)})
    return pd.DataFrame(rows)


# ======================================================================================
def main() -> int:
    pb = pd.read_csv(C.TABLES / "pseudobulk_scores.tsv", sep="\t")
    scored = {
        t: sc.read_h5ad(C.PROC / f"{t}.scored.h5ad")
        for t in ["human_GSE288028", "mouse_GSE288028", "mouse_GSE288529"]
    }

    print("A. zonal gradient", flush=True)
    a_s, a_a = zonal_gradient(pb)
    a_s.to_csv(C.TABLES / "A_zonal_gradient_per_sample.tsv", sep="\t", index=False)
    a_a.to_csv(C.TABLES / "A_zonal_gradient_summary.tsv", sep="\t", index=False)
    print(a_a[a_a.gene_set != "secondary"].to_string(index=False)[:2000], flush=True)

    print("\nB. age effect", flush=True)
    b_tab, gate = age_effect(pb)
    b_tab.to_csv(C.TABLES / "B_age_effect_donor_spread.tsv", sep="\t", index=False)
    (C.TABLES / "B_age_adequacy_gate.json").write_text(json.dumps(gate, indent=2))
    boot = technical_bootstrap_rz(scored)
    boot.to_csv(C.TABLES / "B_resting_zone_TECHNICAL_bootstrap.tsv", sep="\t", index=False)
    print(json.dumps(gate, indent=2), flush=True)

    print("\nC. boundary crossing", flush=True)
    c_tab = boundary_crossing(scored)
    c_tab.to_csv(C.TABLES / "C_boundary_crossing.tsv", sep="\t", index=False)

    print("\nD. GH arm", flush=True)
    hum = scored["human_GSE288028"]
    gp = identify_gp_subclusters(hum)
    gp.to_csv(C.TABLES / "D_rz_subcluster_culture_depletion.tsv", sep="\t", index=False)
    d_all = gh_contrast(hum, restrict_gp2=False)
    d_gp2 = gh_contrast(hum, restrict_gp2=True)
    pd.concat([d_all, d_gp2]).to_csv(C.TABLES / "D_gh_contrast.tsv", sep="\t", index=False)
    hum.write(C.PROC / "human_GSE288028.scored.h5ad", compression="gzip")

    print("\n6. controls", flush=True)
    permutation_control(scored, C.N_PERMUTATIONS).to_csv(
        C.TABLES / "ctrl_permutation_zone_shuffle.tsv", sep="\t", index=False)
    technical_covariate_control(pb).to_csv(
        C.TABLES / "ctrl_technical_covariates.tsv", sep="\t", index=False)
    leave_one_donor_out(pb).to_csv(
        C.TABLES / "ctrl_leave_one_donor_out.tsv", sep="\t", index=False)

    # pathway components by zone
    expr_cols = [c for c in pb.columns if c.startswith("expr_")]
    pw = (
        pb[pb.zone.isin(C.ZONE_ORDER)]
        .groupby(["dataset_tag", "zone"], observed=True)[expr_cols]
        .mean()
        .reset_index()
    )
    pw.to_csv(C.TABLES / "pathway_components_by_zone.tsv", sep="\t", index=False)
    print("\ndone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
