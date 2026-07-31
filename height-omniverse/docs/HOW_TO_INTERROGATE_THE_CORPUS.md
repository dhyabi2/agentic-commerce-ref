# How to interrogate the corpus

```bash
.venv/bin/python -c "import duckdb; con=duckdb.connect('database/height_omniverse.duckdb'); \
  print(con.sql('SHOW TABLES'))"
```
Or open `database/height_omniverse.duckdb` in any DuckDB client. Example queries (all return
source-linked rows):

```sql
-- 1. All interventions with HUMAN adult-height evidence (surrogate actually translated)
SELECT intervention, population, adult_height_effect, source_id
FROM surrogate_calibration WHERE upper(surrogate_translated)='YES';

-- 2. The risky class: positive surrogate but NULL/UNKNOWN final height
SELECT intervention, population, surrogate_effect
FROM surrogate_calibration
WHERE direction='+' AND upper(surrogate_translated) IN ('NO','UNKNOWN');

-- 3. Human-genetics "brake" targets: reduced activity -> taller
SELECT symbol, monogenic_stature_evidence, source_provenance
FROM genes WHERE direction_of_effect ILIKE '%taller%';

-- 4. All evidence edges pointing at ADULT height, with tier
SELECT subject, direction, object, evidence_tier, source_id
FROM evidence_edges WHERE object ILIKE '%adult height%' ORDER BY evidence_tier;

-- 5. Everything recorded for a named compound
SELECT * FROM interventions WHERE name ILIKE '%vosoritide%';
SELECT * FROM surrogate_calibration WHERE intervention ILIKE '%vosoritide%';

-- 6. Contradictory allelic evidence in the same gene (e.g. NPR2 / FGFR3 GoF vs LoF)
SELECT gene, variant_class, height_direction, source_id
FROM read_csv_auto('data/curated/genetic_direction.csv') WHERE gene IN ('NPR2','FGFR3')
ORDER BY gene;

-- 7. Mechanical interventions and whether they measured FINAL length vs a surrogate
SELECT citation, endpoint, direction, directness_tier, source_location
FROM studies WHERE intervention ILIKE ANY (['%PEMF%','%LIPUS%','%compression%','%jump%','%loading%']);

-- 8. Blocked / unretrieved sources (no silent failure)
SELECT * FROM access_failures ORDER BY importance;

-- 9. Genome-wide gene universe size + how many carry recorded evidence
SELECT count(*) AS all_genes,
       count(*) FILTER (WHERE direction_of_effect<>'') AS with_direction FROM genes;
```

DuckDB also queries the Parquet/CSV directly without the DB, e.g.
`SELECT * FROM 'data/processed/ALL_GENES.parquet' LIMIT 5;`

Knowledge graph: `knowledge_graph/height_omniverse.graphml` (open in Cytoscape/Gephi) or
`nodes.parquet`/`edges.parquet` for Neo4j import. Dashboard: `.venv/bin/streamlit run dashboard/app.py`.
