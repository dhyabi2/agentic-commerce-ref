# Deviations from the pre-registration

Every departure from `docs/PREREGISTRATION.md`, with its reason.

## 1. Non-chondrocyte clusters excluded from zonal analysis

**Deviation.** The pre-registration assigns every cluster to one of four zones. The
implementation adds a fifth label, `non_chondrocyte`.

**Reason.** The growth plate digests contain immune, endothelial and erythroid cells.
Forcing those clusters into a chondrocyte zone would corrupt every zonal statistic.
Clusters are labelled `non_chondrocyte` when their mean chondrocyte score
(COL2A1/ACAN/SOX9/COL9A1/COL11A1/COL11A2) is <= 0, or when their immune (PTPRC/CD68/LYZ/CD3E),
endothelial (PECAM1/CDH5/KDR) or erythroid (HBB/HBA1/ALAS2) score exceeds it. These cells
are excluded from all zonal and pseudobulk inference and are reported in the cell-count
tables. The gate uses lineage markers only — none of the scoring genes.

## 2. Zone and activity scoring use the unregressed expression layer

**Deviation.** The brief instructs that cell cycle and stress signatures be regressed out.
Regression is applied to the highly-variable-gene matrix used for PCA, neighbours,
clustering and UMAP. Zone-marker scoring and HDAC4 activity scoring are computed on the
**unregressed** log-normalised layer.

**Reason.** The proliferative zone is *defined* by CCND1/MKI67/E2F1. Scoring those markers
on a matrix from which the cell-cycle signature has been removed would erase the very
signal used to identify the zone. Regressing the embedding while scoring the unregressed
values keeps the instruction's intent (cell cycle should not drive cluster structure)
without making the proliferative zone unrecoverable. Both matrices are retained in the
`.h5ad` (`layers["counts"]`, `layers["lognorm"]`) so either choice can be re-derived.

## 3. `regress_out` implemented as a single least-squares solve

**Deviation.** `scanpy.pp.regress_out` replaced by `regress_out_fast` in `src/02_process.py`.

**Reason.** Performance only. Scanpy loops per gene; the replacement solves the same OLS
problem for all genes at once via `numpy.linalg.lstsq` against the same design matrix
(intercept + S_score + G2M_score + stress_score + pct_counts_mt). The residuals are
mathematically identical.

## 4. Supplementary rat dataset added (analysis E)

**Deviation.** GSE16981 (rat, bulk Affymetrix, microdissected zones, ages 1/3/6/9/12 wk)
was not named in the brief, which asked for growth plate **scRNA-seq** with an age series.

**Reason.** The GEO search (section 2 of the brief, `results/tables/geo_search_log.tsv`)
established that **no** growth plate scRNA-seq dataset with a developmental age series
exists in GEO. GSE16981 is the only public growth plate dataset with genuine zone x age
structure and biological replication, so it is the only available empirical test of the
temporal prediction. It is reported as clearly-labelled supplementary analysis E, is not
part of the primary result, and carries its own limitations (rat; bulk; the 3-12 wk
animals were castrated at day 18).

## 5. Analysis B not fitted

**Not a deviation — this is the pre-registered adequacy gate firing.**

Gate (a) of PREREGISTRATION section 5 requires per-donor developmental age to exist in the
metadata. It does not, for any dataset analysed. As pre-registered, no age model was
fitted and no surrogate regressor was substituted. What is reported instead is the
observed between-donor spread in the resting zone, which bounds what any age effect would
have had to exceed to be detectable at n=4.
