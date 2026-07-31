#!/usr/bin/env python
"""SUPPLEMENTARY analysis E: the only public growth plate dataset with zone x age structure.

GSE16981 (Lui/Andrade/Nilsson/Baron) -- rat proximal tibial growth plate, zones separated
by manual microdissection, Affymetrix Rat 230 2.0, MAS5 signal, n=5 biological replicates
per group, ages 1, 3, 6, 9 and 12 weeks.

This is NOT the primary analysis and does NOT substitute for it:
  * rat, not human -- and rat growth plates do not undergo the human fusion event
  * bulk microarray, not scRNA-seq -- no single-cell resolution, no zone inference
  * animals older than 1 wk were CASTRATED at 18 days, so the 1wk-vs-older contrast is
    confounded with loss of gonadal steroids. The 3/6/9/12wk series is internally
    consistent (all castrated) and is the contrast used for the age model.
  * transcript abundance still cannot measure HDAC4 subcellular localisation.

Outputs: results/tables/E_rat_*.tsv
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hdac4_config as C  # noqa: E402

MATRIX = C.RAW / "GSE16981_series_matrix.txt.gz"
ANNOT = C.RAW / "GPL1355.annot.gz"

ZONE_CANON = {
    "resting zone": "resting",
    "proliferative zone": "proliferative",
    "hypertrophic zone": "hypertrophic",
    "early hypertrophic zone": "prehypertrophic",
}


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_lines, data_lines, in_table = [], [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            (data_lines if in_table else meta_lines).append(line)

    gsms = None
    chars: list[list[str]] = []
    for line in meta_lines:
        if line.startswith("!Sample_geo_accession"):
            gsms = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip('"') for x in line.split("\t")[1:]])

    obs = pd.DataFrame({"gsm": gsms})
    for row in chars:
        keys = {c.split(":", 1)[0].strip().lower() for c in row if ":" in c}
        if len(keys) == 1:
            k = keys.pop()
            obs[k] = [c.split(":", 1)[1].strip() if ":" in c else np.nan for c in row]

    hdr = [x.strip('"') for x in data_lines[0].split("\t")]
    recs = [l.split("\t") for l in data_lines[1:] if l.strip()]
    expr = pd.DataFrame(recs, columns=hdr).set_index(hdr[0])
    expr.index = [i.strip('"') for i in expr.index]
    expr = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    expr.columns = [c.strip('"') for c in expr.columns]
    return obs, expr


def probe_to_symbol(path: Path) -> pd.Series:
    rows, started = [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                started = True
                header = next(fh).rstrip("\n").split("\t")
                i_id, i_sym = header.index("ID"), header.index("Gene symbol")
                continue
            if line.startswith("!platform_table_end"):
                break
            if started:
                p = line.rstrip("\n").split("\t")
                if len(p) > max(i_id, i_sym) and p[i_sym]:
                    rows.append((p[i_id], p[i_sym].split("///")[0].strip()))
    return pd.Series(dict(rows))


def collapse_to_genes(expr: pd.DataFrame, mapping: pd.Series) -> pd.DataFrame:
    """log2 MAS5, median-centre arrays, then take the max-expressed probe per gene."""
    x = np.log2(expr.clip(lower=1.0))
    x = x - x.median(axis=0) + x.median(axis=0).mean()  # per-array median scaling
    sym = mapping.reindex(x.index)
    x = x[sym.notna().to_numpy()]
    sym = sym.dropna()
    x = x.assign(_sym=sym.values, _mean=x.mean(1).values)
    x = x.sort_values("_mean", ascending=False).drop_duplicates("_sym")
    x = x.set_index("_sym").drop(columns=["_mean"])
    x.index.name = "gene"
    return x


def zscore_score(gexpr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    """Mean of per-gene z-scores across samples -- the bulk analogue of score_genes."""
    found = [g for g in genes if g in gexpr.index]
    if not found:
        return pd.Series(np.nan, index=gexpr.columns), []
    sub = gexpr.loc[found]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1, ddof=0).replace(0, np.nan), axis=0)
    return z.mean(0), found


def ulm_score(gexpr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """decoupler ULM on the bulk matrix (samples x genes)."""
    import decoupler as dc

    found = [g for g in genes if g in gexpr.index]
    if not found:
        return pd.Series(np.nan, index=gexpr.columns)
    mat = gexpr.T  # samples x genes
    net = pd.DataFrame({"source": "MEF2C_RUNX2", "target": found, "weight": 1.0})
    try:
        res = dc.mt.ulm(mat, net, tmin=min(3, len(found)), verbose=False)
        est = res[0] if isinstance(res, tuple) else res
        return pd.Series(np.asarray(est).ravel()[: mat.shape[0]], index=mat.index)
    except Exception as exc:  # noqa: BLE001
        print(f"  [rat ULM failed: {exc}]")
        return pd.Series(np.nan, index=gexpr.columns)


def main() -> int:
    obs, expr = parse_series_matrix(MATRIX)
    mapping = probe_to_symbol(ANNOT)
    gexpr = collapse_to_genes(expr, mapping)
    print(f"parsed {expr.shape[0]} probes -> {gexpr.shape[0]} genes x {gexpr.shape[1]} arrays")

    obs["age_weeks"] = obs["age"].str.extract(r"(\d+)").astype(float)
    obs["zone"] = (
        obs["growth plate zone"].str.strip().str.lower().map(ZONE_CANON)
        .fillna(obs["growth plate zone"].str.strip().str.lower())
    )
    obs["castrated"] = obs.get("castrated", pd.Series("unknown", index=obs.index)).str.lower()
    obs = obs.set_index("gsm").loc[gexpr.columns]
    print(obs.groupby(["age_weeks", "zone", "castrated"], observed=True).size().to_string())

    # scores (rat symbols follow the same capitalisation rule as mouse)
    rat = lambda gs: C.map_genes(gs, "mouse")  # noqa: E731
    for name, genes in [
        ("primary", C.TARGETS_PRIMARY),
        ("secondary", C.TARGETS_SECONDARY),
        ("housekeeping", C.HOUSEKEEPING),
    ]:
        s, found = zscore_score(gexpr, rat(genes))
        obs[f"target_{name}_zmean"] = s.reindex(obs.index).values
        obs[f"target_{name}_ulm"] = ulm_score(gexpr, rat(genes)).reindex(obs.index).values
        obs[f"hdac4_activity_{name}_zmean"] = -obs[f"target_{name}_zmean"]
        print(f"  {name}: {len(found)}/{len(genes)} genes found -> {found}")

    for g in rat(C.PATHWAY_COMPONENTS):
        if g in gexpr.index:
            obs[f"expr_{g}"] = gexpr.loc[g].reindex(obs.index).values

    obs.reset_index().to_csv(C.TABLES / "E_rat_sample_scores.tsv", sep="\t", index=False)

    # ---- E1: zonal gradient at 1 wk (all three zones present, uncastrated) ----
    rows = []
    one = obs[obs.age_weeks == 1]
    rank = {"resting": 0, "proliferative": 1, "prehypertrophic": 2, "hypertrophic": 3}
    for col in ["target_primary_zmean", "target_secondary_zmean", "target_housekeeping_zmean",
                "target_primary_ulm", "target_secondary_ulm", "target_housekeeping_ulm"]:
        d = one[one.zone.isin(rank)].dropna(subset=[col])
        if d.zone.nunique() < 3:
            continue
        r = stats.spearmanr(d.zone.map(rank), d[col])
        means = d.groupby("zone", observed=True)[col].mean().to_dict()
        rows.append({"score": col, "n_arrays": len(d), "spearman_rho": r.statistic,
                     "p_value": r.pvalue,
                     **{f"mean_{k}": round(float(v), 4) for k, v in means.items()}})
    pd.DataFrame(rows).to_csv(C.TABLES / "E_rat_zonal_gradient_1wk.tsv", sep="\t", index=False)

    # ---- E2: age effect within zone, castrated animals only (3/6/9/12 wk) ----
    rows = []
    for zone in sorted(obs.zone.dropna().unique()):
        for col in ["target_primary_zmean", "target_secondary_zmean",
                    "target_housekeeping_zmean", "target_primary_ulm"]:
            # The 1 wk arm cannot be pooled into the age series: those animals were
            # dissected into SEPARATE resting/proliferative/hypertrophic zones, whereas the
            # 3-12 wk animals give a single combined "proliferative and early hypertrophic"
            # fraction -- and the 1 wk animals were not castrated. Different tissue
            # compartment and different gonadal status, so it is reported separately.
            for label, sel in [
                ("castrated_only", obs[(obs.zone == zone) & (obs.castrated == "yes")]),
                ("all_castration_states", obs[obs.zone == zone]),
            ]:
                d = sel.dropna(subset=[col, "age_weeks"])
                if d.age_weeks.nunique() < 3 or len(d) < 6:
                    continue
                x, y = d.age_weeks.to_numpy(float), d[col].to_numpy(float)
                lr = stats.linregress(x, y)
                n = len(x)
                tcrit = stats.t.ppf(0.975, n - 2)
                sp = stats.spearmanr(x, y)
                rows.append({
                    "zone": zone, "score": col, "subset": label,
                    "n_arrays": n, "n_ages": int(d.age_weeks.nunique()),
                    "ages": ",".join(str(int(a)) for a in sorted(d.age_weeks.unique())),
                    "slope_per_week": lr.slope, "stderr": lr.stderr,
                    "ci95_lo": lr.slope - tcrit * lr.stderr,
                    "ci95_hi": lr.slope + tcrit * lr.stderr,
                    "p_value": lr.pvalue, "r_squared": lr.rvalue ** 2,
                    "spearman_rho": sp.statistic, "spearman_p": sp.pvalue,
                    "all_castrated": bool((d.castrated == "yes").all()),
                    "zone_dissection_consistent": bool(d.zone.nunique() == 1),
                    "interpretation_if_significant": (
                        "target expression RISES with age -> inferred HDAC4 repressive "
                        "activity FALLS with age (direction predicted by the hypothesis)"
                        if lr.slope > 0 else
                        "target expression FALLS with age -> inferred HDAC4 repressive "
                        "activity RISES with age (opposite to the prediction)"
                    ),
                })
    pd.DataFrame(rows).to_csv(C.TABLES / "E_rat_age_effect.tsv", sep="\t", index=False)

    # ---- E3b: leave-one-age-out on the age slope (is it carried by one age group?) ----
    rows = []
    for zone in sorted(obs.zone.dropna().unique()):
        base = obs[(obs.zone == zone) & (obs.castrated == "yes")]
        if base.age_weeks.nunique() < 4:
            continue
        for col in ["target_primary_zmean", "target_secondary_zmean",
                    "target_housekeeping_zmean", "target_primary_ulm"]:
            for drop in [None, *sorted(base.age_weeks.unique())]:
                d = base if drop is None else base[base.age_weeks != drop]
                d = d.dropna(subset=[col])
                if d.age_weeks.nunique() < 3:
                    continue
                lr = stats.linregress(d.age_weeks.to_numpy(float), d[col].to_numpy(float))
                rows.append({
                    "zone": zone, "score": col,
                    "age_left_out_weeks": "none (full)" if drop is None else int(drop),
                    "n_arrays": len(d), "slope_per_week": lr.slope,
                    "p_value": lr.pvalue, "r_squared": lr.rvalue ** 2,
                })
    pd.DataFrame(rows).to_csv(C.TABLES / "E_rat_leave_one_age_out.tsv", sep="\t", index=False)

    # ---- E3: per age x zone means ----
    score_cols = [c for c in obs.columns if c.startswith(("target_", "hdac4_", "expr_"))]
    (obs.groupby(["zone", "age_weeks", "castrated"], observed=True)[score_cols]
        .agg(["mean", "std", "size"])
        .to_csv(C.TABLES / "E_rat_zone_age_means.tsv", sep="\t"))

    print("\ndone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
