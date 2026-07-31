# Pth1r across growth-plate chondrocyte states — GSE288529

**Question.** Is *Pth1r* expression low in *Apoe⁺/Clu⁺/C1qtnf3⁺* resting-zone chondrocytes and
higher in *Ihh⁺/Runx2⁺* prehypertrophic chondrocytes?

**Short answer.** Yes. Independently of *Pth1r*, resting and prehypertrophic states were defined
from their canonical markers; *Pth1r* is then **~5× higher per cell** and **~15× higher in
pseudobulk** in prehypertrophic than in resting chondrocytes, with the effect visible across every
metric (mean, median, % detected, pseudobulk) and overwhelmingly significant. The pattern is exactly
the textbook PTHrP→PTH1R (Ihh–PTHrP) feedback loop, which is a strong internal-validity check.

---

## Data & methods

- **Dataset:** GSE288529, `GSE288529_RAW.tar` from NCBI GEO (downloaded programmatically).
  Single 10x sample **GSM8769462** (mouse, `filtered_feature_bc_matrix.h5`), 12,636 called cells ×
  32,286 features, raw UMI counts.
- **Environment:** isolated Python venv (`.venv/`), Python 3.11.15, Scanpy 1.11.5, numpy 2.4.6,
  scipy 1.17.1, pandas 2.3.3, leidenalg. System Python untouched. Full pins in
  `environment_requirements.txt`.
- **Pipeline** (`scripts/01_process.py`): QC (min 200 genes/cell, genes in ≥3 cells, mito < 15%,
  drop top 0.5% total-count outliers → **11,576 cells × 24,697 genes**); preserve raw counts in a
  `counts` layer; `normalize_total(1e4)` + `log1p`; 2,000 HVGs; scale; PCA(50); neighbors(15,
  30 PCs); UMAP; Leiden (res 1.0) → 24 clusters.
- **Annotation** (`scripts/02_assign.py`, `scripts/03_annotate_quantify.py`): each cluster scored
  against marker panels using the mean of per-gene z-scores across clusters; **Pth1r was never used
  for clustering or annotation.** Data-driven scores were then reconciled with the raw per-cluster
  marker table and biological judgment; small low-count clusters with only noise-level scores were
  set aside rather than forced into a state.

### Cluster → state map (Pth1r-independent)

| State | Leiden clusters | Defining evidence |
|---|---|---|
| **Resting** | 8 | Highest *Apoe* (3.2, 97% det) + *Clu* (1.3, 77%); lowest hypertrophic markers |
| **Proliferative** | 4, 5, 16 | High *Acan/Col2a1/Sox9*, highest *Ccnd1*, *C1qtnf3⁺*, *Ihh* low |
| **Prehypertrophic** | 3, 7 | *Ihh⁺* (56–91% det) with *Runx2* (97%) + *Sp7* (82–86%); *Mmp13* low |
| **Hypertrophic** | 6, 12, 17, 19 | High *Col10a1* and/or *Mmp13/Ibsp* (terminal) |
| *Excluded — Articular* | 9, 10, 15 | *Prg4*(up to 100% det)/*Cytl1* high |
| *Excluded — Osteo/Perichondrial* | 0, 14 | *Col1a1/Bglap/Ibsp* (0) or *Col1a1/Col3a1/Postn* (14) |
| *Excluded — Endothelial* | 1 | *Cdh5/Pecam1/Kdr* (>90% det) |
| *Excluded — Immune* | 13 | *Ptprc⁺* |
| *Excluded — LowQual/Mixed/border* | 2, 11, 18, 20, 21, 22, 23 | No dominant panel; mixed/low-count |

Cluster 6 (*Ihh⁺, Col10a1* ~99%) is a genuine prehypertrophic→hypertrophic transition; it was placed
in Hypertrophic on *Col10a1* dominance. Reassigning it either way does not change the conclusion.

---

## Result 1 — Pth1r per state (`tables/pth1r_by_state.csv`)

| State | n cells | mean lognorm | median lognorm | % detected | pseudobulk (raw/10k) |
|---|---:|---:|---:|---:|---:|
| Resting | 914 | **0.550** | 0.484 | **63.3%** | **1.15** |
| Proliferative | 2,963 | 0.918 | 0.845 | 83.8% | 2.80 |
| **Prehypertrophic** | 1,576 | **2.683** | 2.862 | **98.7%** | **17.0** |
| Hypertrophic | 1,624 | 1.491 | 1.419 | 62.9% | 16.2 |

*Pth1r* rises monotonically resting → proliferative → prehypertrophic and peaks at the prehypertrophic
stage on a per-cell basis. (Hypertrophic cells have high pseudobulk because they are large/high-count,
but a lower detection fraction — consistent with receptor being maximal at the prehypertrophic exit.)

## Result 2 — Resting vs Prehypertrophic (`tables/resting_vs_prehyp.csv`)

| Metric | Value | 95% CI (bootstrap, 5,000×) |
|---|---:|---|
| Mean lognorm (Resting / Prehyp) | 0.550 / 2.683 | — |
| **Expression ratio (Pre/Rest)** | **4.88×** | 4.56 – 5.22 |
| Mean-difference (Pre − Rest) | 2.13 | 2.08 – 2.18 |
| % detected (Resting / Prehyp) | 63.3% / 98.7% | — |
| **Detection ratio (Pre/Rest)** | **1.56×** | 1.49 – 1.64 |
| **Pseudobulk ratio (Pre/Rest)** | **14.7×** | — |
| Wilcoxon rank-sum | z = 39.5, **p < 1e-300** | — |
| Mann–Whitney U | U = 1.40e6, **p < 1e-300** | — |

Every metric agrees and every confidence interval excludes the null (ratio = 1 / difference = 0).

## Result 3 — Are Pth1r⁺ resting cells a distinct subpopulation? (`tables/resting_subpopulation.csv`)

Tested on the resting-only kNN graph (Pth1r never used to build it):

- *Pth1r⁺* resting fraction: **63.3%**.
- **kNN homophily** of *Pth1r⁺* cells: observed **0.671** vs permutation null 0.633 ± 0.008
  → **z = 4.76, p ≈ 5e-4**.
- **Moran's I** of *Pth1r* within resting: **0.185** vs null ≈ 0 → **p ≈ 5e-4**.

**Interpretation:** *Pth1r⁺* resting cells are **not randomly scattered** — there is statistically
significant, but *modest*, spatial structure (Moran's I ≈ 0.19). Practically, this is a graded
sub-signal: the resting cells carrying more *Pth1r* sit preferentially toward the resting↔proliferative
border rather than forming a sharply isolated island. So *Pth1r* is broadly (if weakly) present across
the reserve pool, with a mild gradient toward the transit-amplifying edge.

---

## scRNA-seq dropout — does it explain the gap, and is the claim supported?

**Dropout (technical zero-inflation) is real here and it matters most for the low-expressing resting
group.** Three points:

1. **Direction of bias.** Dropout inflates zeros preferentially where transcripts are scarce, i.e. in
   *resting* cells. So the *measured* resting level is, if anything, a slight **under-estimate** of the
   true value. This makes the resting–prehypertrophic gap a conservative floor, not an artifact.

2. **Concordance of independent estimators.** If the difference were a dropout artifact it would show
   up in the detection fraction but shrink in count-based measures. The opposite is seen: detection
   ratio is the *smallest* effect (1.56×), per-cell mean is larger (4.9×), and **library-size-normalized
   pseudobulk is largest (14.7×)**. Count-based, dropout-robust pseudobulk *amplifies* the gap. That
   ordering is what a real biological difference looks like, not a technical one.

3. **What dropout does compress.** Because 63% of resting cells still show ≥1 *Pth1r* UMI, resting
   chondrocytes are **not receptor-null** — they express the receptor at a genuinely lower level, not
   zero. The per-cell mean ratio (~5×) is partly a mixture of "fewer positive cells" and "lower level
   when positive," and dropout blurs the exact per-cell magnitude. The robust, defensible statement is
   the *ranking and the pseudobulk fold-change*, not a precise single-cell copy number.

### Verdict on the claim

> "Systemic PTH1R agonism primarily targets the prehypertrophic exit rather than the resting-zone
> reserve because Pth1r expression is much lower in resting cells."

**The expression premise is well supported; the pharmacological conclusion is a reasonable but
partial extrapolation.**

- **Supported:** *Pth1r* mRNA is unambiguously and substantially lower in resting than in
  prehypertrophic chondrocytes — 5× lower per cell, ~15× lower in pseudobulk, p < 1e-300, robust to
  dropout, and matching the canonical Ihh–PTHrP loop (PTHrP source at the resting/periarticular end,
  PTH1R receptor peak at the prehypertrophic exit). A systemic agonist will therefore encounter far
  more receptor-bearing target in the prehypertrophic zone.
- **Caveats before treating it as settled:**
  1. This is **one animal / one 10x library** — no biological replication, so no cross-sample CI on
     the fold-change. The within-sample statistics are strong but n(samples)=1.
  2. **mRNA ≠ receptor protein ≠ signaling output.** Agonist effect depends on surface receptor and
     downstream cAMP/PKA, which are not measured here.
  3. Resting cells are **not receptor-negative** (63% detected). "Primarily targets prehypertrophic"
     is directionally right, but the reserve zone is not pharmacologically silent — a lower-magnitude
     resting response is plausible.
  4. "Exit" is inferred from a static snapshot; no lineage/pseudotime causality is established here.

**Bottom line:** the data strongly back *"Pth1r expression is much lower in resting cells"* and are
consistent with prehypertrophic cells being the dominant transcriptional target of PTH1R agonism.
They do **not**, on their own, prove the functional/pharmacological claim — that needs protein-level
readouts and biological replicates.

---

## Files

```
growth_plate_pth1r/
├── REPORT.md                       # this file
├── environment_requirements.txt    # exact package pins
├── scripts/
│   ├── 01_process.py               # QC → normalize → PCA → UMAP → Leiden
│   ├── 02_assign.py                # z-scored panel assignment (Pth1r-free)
│   └── 03_annotate_quantify.py     # annotate, quantify, stats, subpop test, figures
├── tables/
│   ├── pth1r_by_state.csv          # main summary table (Result 1)
│   ├── resting_vs_prehyp.csv       # Result 2 (ratios, CIs, tests)
│   ├── resting_subpopulation.csv   # Result 3 (homophily, Moran's I)
│   ├── cluster_marker_means.csv    # per-cluster mean expr (annotation evidence)
│   ├── cluster_marker_pctdetected.csv
│   └── cluster_assignment_scores.csv
├── figures/
│   ├── umap_states.png, umap_leiden.png
│   ├── umap_Apoe.png, umap_Pthlh.png, umap_Ihh.png, umap_Col10a1.png, umap_Pth1r.png
│   ├── dotplot_markers_by_cluster.png, dotplot_Pth1r_by_state.png
│   └── violin_Pth1r_by_state.png
└── data/  (GSE288529_RAW.tar, extracted 10x files, processed.h5ad)  # gitignored if large
```

### Reproduce
```bash
python3 -m venv .venv && ./.venv/bin/pip install -r environment_requirements.txt
# place GSE288529_RAW.tar in data/ and extract, then:
./.venv/bin/python scripts/01_process.py
./.venv/bin/python scripts/02_assign.py
./.venv/bin/python scripts/03_annotate_quantify.py
```
