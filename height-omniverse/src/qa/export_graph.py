#!/usr/bin/env python
"""Export the knowledge graph (evidence_edges) to GraphML + nodes/edges Parquet (Neo4j-friendly)."""
import duckdb, pandas as pd, pathlib, networkx as nx
ROOT = pathlib.Path(__file__).resolve().parents[2]
con = duckdb.connect(str(ROOT / "database" / "height_omniverse.duckdb"))
edges = con.execute("SELECT * FROM evidence_edges").fetchdf()
con.close()
KG = ROOT / "knowledge_graph"; KG.mkdir(exist_ok=True)
if edges.empty:
    print("no edges yet"); raise SystemExit(0)
G = nx.MultiDiGraph()
for _, e in edges.iterrows():
    for n in (e["subject"], e["object"]):
        if not G.has_node(n): G.add_node(n, label=str(n))
    G.add_edge(e["subject"], e["object"], key=e["edge_id"],
               predicate=str(e["predicate"]), direction=str(e["direction"]),
               context=str(e["context"]), tier=str(e["evidence_tier"]),
               source_id=str(e["source_id"]), confidence=str(e["confidence"]))
nx.write_graphml(G, KG / "height_omniverse.graphml")
nodes = pd.DataFrame({"node_id": list(G.nodes())})
nodes.to_parquet(KG / "nodes.parquet", index=False)
edges.to_parquet(KG / "edges.parquet", index=False)
print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges -> {KG}")
