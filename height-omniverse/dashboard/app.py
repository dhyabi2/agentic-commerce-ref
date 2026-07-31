#!/usr/bin/env python
"""Minimal Streamlit explorer over the Height Omniverse DuckDB warehouse.
Run: .venv/bin/streamlit run dashboard/app.py
Clearly separates EMPIRICAL evidence, INFERENCE, and SIMULATION. Every row shows a source_id."""
import duckdb, pathlib, streamlit as st, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[1]
con = duckdb.connect(str(ROOT / "database" / "height_omniverse.duckdb"), read_only=True)

st.set_page_config(page_title="Height Omniverse", layout="wide")
st.title("Height Omniverse — evidence warehouse")
st.warning("Evidence & simulation warehouse. NOT a treatment recommendation or personal protocol. "
           "No dosing/timing/'best stack' output.")

tab = st.sidebar.radio("View", ["Genes", "Interventions", "Surrogate calibration",
                                "Evidence edges", "Studies", "Simulation (twin)"])
def show(df): st.dataframe(df, use_container_width=True, height=560)

if tab == "Genes":
    q = st.text_input("Gene symbol / alias contains")
    sql = "SELECT symbol,name,locus_group,direction_of_effect,monogenic_stature_evidence,source_provenance FROM genes"
    if q: sql += f" WHERE symbol ILIKE '%{q}%' OR aliases ILIKE '%{q}%' OR monogenic_stature_evidence <> ''"
    else: sql += " WHERE monogenic_stature_evidence <> '' OR direction_of_effect <> ''"
    show(con.execute(sql + " LIMIT 500").fetchdf())
    st.caption(f"Gene universe size: {con.execute('SELECT count(*) FROM genes').fetchone()[0]}")
elif tab == "Interventions":
    show(con.execute("SELECT name,modality,target,action,approval_status,route,pk_summary FROM interventions").fetchdf())
elif tab == "Surrogate calibration":
    st.subheader("Did the surrogate translate to ATTAINED adult height?")
    df = con.execute("SELECT * FROM surrogate_calibration").fetchdf()
    show(df)
    st.caption("YES = surrogate banked as adult height; NO = surrogate gain did NOT translate; "
               "UNKNOWN = final-height data absent.")
elif tab == "Evidence edges":
    show(con.execute("SELECT subject,predicate,object,direction,context,evidence_tier,confidence,source_id FROM evidence_edges").fetchdf())
elif tab == "Studies":
    show(con.execute("SELECT citation,species,age_stage,disease,intervention,endpoint,effect_size,direction,directness_tier FROM studies").fetchdf())
elif tab == "Simulation (twin)":
    st.info("SIMULATION output — assumption-exposing, NOT evidence.")
    p = ROOT / "results" / "models" / "growth_deltaL_by_family_stage.csv"
    if p.exists(): show(pd.read_csv(p))
    fig = ROOT / "results" / "figures" / "deltaL_posteriors.png"
    if fig.exists(): st.image(str(fig))
