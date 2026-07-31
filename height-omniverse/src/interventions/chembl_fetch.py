#!/usr/bin/env python
"""Phase F (reproducible): fetch compound-target bioactivity from the ChEMBL REST API for height
targets, writing results/tables/COMPOUND_TARGET_full.parquet. Run to expand/refresh the layer.
No efficacy is inferred from binding. Requires network."""
import requests, pandas as pd, pathlib, time
ROOT = pathlib.Path(__file__).resolve().parents[2]
BASE = "https://www.ebi.ac.uk/chembl/api/data"
TARGETS = {"GHR":"growth hormone receptor","IGF1R":"IGF1R","NPR2":"NPR2","FGFR3":"FGFR3","PTH1R":"PTH1R",
           "ESR1":"estrogen receptor alpha","CYP19A1":"aromatase","AR":"androgen receptor","IGF1":"IGF1",
           "PAPPA2":"pappalysin 2","PDE3A":"PDE3A","SMO":"smoothened","CTSK":"cathepsin K","MMP13":"MMP13"}
def get(url, params):
    r = requests.get(url, params=params, headers={"Accept":"application/json"}, timeout=60); r.raise_for_status(); return r.json()
rows=[]
for gene,q in TARGETS.items():
    try:
        ts=get(f"{BASE}/target/search",{"q":q,"format":"json"}).get("targets",[])
        tid=next((t["target_chembl_id"] for t in ts if t.get("organism")=="Homo sapiens" and t.get("target_type")=="SINGLE PROTEIN"),None)
        if not tid: continue
        acts=get(f"{BASE}/activity",{"target_chembl_id":tid,"standard_type__in":"IC50,Ki,Kd,EC50","limit":200,"format":"json"}).get("activities",[])
        for a in acts:
            if a.get("standard_value") and a.get("standard_units")=="nM":
                rows.append({"target_gene":gene,"target_chembl_id":tid,"molecule_chembl_id":a.get("molecule_chembl_id"),
                             "standard_type":a.get("standard_type"),"standard_value_nM":a.get("standard_value"),
                             "assay_chembl_id":a.get("assay_chembl_id"),"source":"ChEMBL API"})
        time.sleep(0.3)
    except Exception as e:
        print("skip",gene,e)
pd.DataFrame(rows).to_parquet(ROOT/"results/tables/COMPOUND_TARGET_full.parquet",index=False)
print(f"fetched {len(rows)} activity rows")
