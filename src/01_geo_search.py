#!/usr/bin/env python
"""Search GEO for growth plate scRNA-seq datasets carrying a developmental age series.

Logs every query and its hit count to results/tables/geo_search_log.tsv, and every
distinct series hit to results/tables/geo_search_hits.tsv.

The age range available across public datasets is the binding constraint on the whole
analysis, so this step is a first-class result, not a preliminary.
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Query terms specified in the brief, plus scRNA-seq-restricted combinations.
QUERIES = [
    'growth plate[All Fields]',
    'physis[All Fields]',
    'epiphyseal[All Fields]',
    'resting zone[All Fields]',
    'chondrocyte[All Fields]',
    'growth plate[All Fields] AND "single cell"[All Fields]',
    'growth plate[All Fields] AND (scRNA-seq[All Fields] OR "single cell RNA"[All Fields])',
    'physis[All Fields] AND "single cell"[All Fields]',
    'epiphyseal[All Fields] AND "single cell"[All Fields]',
    'resting zone[All Fields] AND "single cell"[All Fields]',
    'chondrocyte[All Fields] AND "single cell"[All Fields]',
    'growth plate[All Fields] AND postnatal[All Fields]',
    'growth plate[All Fields] AND juvenile[All Fields]',
    'growth plate[All Fields] AND pubertal[All Fields]',
    'growth plate[All Fields] AND (aging[All Fields] OR ageing[All Fields] OR "age series"[All Fields])',
    'chondrocyte[All Fields] AND "single cell"[All Fields] AND (postnatal[All Fields] OR juvenile[All Fields] OR pubertal[All Fields])',
    'growth plate[All Fields] AND "single cell"[All Fields] AND Homo sapiens[Organism]',
    'growth plate[All Fields] AND "single cell"[All Fields] AND Mus musculus[Organism]',
]

# Keywords suggesting the series spans more than one developmental stage.
AGE_SERIES_HINT = re.compile(
    r"\b(P\d{1,3}\b|E\d{1,2}\.?\d?\b|\d+\s*-\s*\d+\s*(week|wk|month|mo|year|yr)|"
    r"week[s]?\s*old|month[s]?\s*old|year[s]?\s*old|postnatal day|"
    r"time course|timecourse|developmental stage|aged?\s+\d|"
    r"juvenile|adolescent|pubertal|adult|neonat|embryonic)\b",
    re.IGNORECASE,
)


def eutils_get(endpoint: str, params: dict, retries: int = 4) -> str:
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except Exception as exc:  # network flake -> exponential backoff
            last = exc
            time.sleep(2 ** (i + 1))
    raise RuntimeError(f"eutils failed after {retries} attempts: {url}: {last}")


def esearch(term: str, retmax: int = 200) -> tuple[int, list[str]]:
    txt = eutils_get(
        "esearch.fcgi",
        {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"},
    )
    js = json.loads(txt)["esearchresult"]
    return int(js.get("count", 0)), js.get("idlist", [])


def esummary(uids: list[str]) -> list[dict]:
    out = []
    for i in range(0, len(uids), 100):
        chunk = uids[i : i + 100]
        txt = eutils_get(
            "esummary.fcgi",
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"},
        )
        js = json.loads(txt).get("result", {})
        for uid in js.get("uids", []):
            out.append(js[uid])
        time.sleep(0.4)
    return out


def main() -> int:
    log_rows, hit_rows, seen = [], [], {}

    for term in QUERIES:
        try:
            count, uids = esearch(term)
        except Exception as exc:
            log_rows.append({"query": term, "n_hits": "ERROR", "note": str(exc)[:200]})
            print(f"[ERROR] {term}: {exc}", file=sys.stderr)
            continue
        log_rows.append({"query": term, "n_hits": count, "note": f"retrieved {len(uids)} uids"})
        print(f"{count:6d}  {term}")
        for rec in esummary(uids):
            acc = rec.get("accession", "")
            if not acc.startswith("GSE") or acc in seen:
                continue
            seen[acc] = True
            title = rec.get("title", "")
            summary = rec.get("summary", "")
            blob = f"{title} {summary}"
            hit_rows.append(
                {
                    "accession": acc,
                    "organism": rec.get("taxon", ""),
                    "n_samples": rec.get("n_samples", ""),
                    "gdstype": rec.get("gdstype", ""),
                    "pdat": rec.get("pdat", ""),
                    "title": title,
                    "is_single_cell": bool(
                        re.search(r"single[- ]cell|scRNA|snRNA|10x|10X", blob, re.I)
                    ),
                    "growth_plate_relevant": bool(
                        re.search(
                            r"growth plate|physis|epiphys|resting zone|chondrocyte", blob, re.I
                        )
                    ),
                    "age_series_hint": bool(AGE_SERIES_HINT.search(blob)),
                    "summary": summary.replace("\t", " ").replace("\n", " ")[:1200],
                }
            )
        time.sleep(0.4)

    pd.DataFrame(log_rows).to_csv(TABLES / "geo_search_log.tsv", sep="\t", index=False)
    hits = pd.DataFrame(hit_rows)
    if not hits.empty:
        hits = hits.sort_values(
            ["growth_plate_relevant", "is_single_cell", "age_series_hint"], ascending=False
        )
    hits.to_csv(TABLES / "geo_search_hits.tsv", sep="\t", index=False)

    cand = hits[
        hits.growth_plate_relevant & hits.is_single_cell
    ] if not hits.empty else hits
    cand.to_csv(TABLES / "geo_search_candidates.tsv", sep="\t", index=False)

    print(f"\nqueries run          : {len(QUERIES)}")
    print(f"distinct GSE hits    : {len(hits)}")
    print(f"growth-plate + scRNA : {len(cand)}")
    if not cand.empty:
        print(f"  ...with age hint   : {int(cand.age_series_hint.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
