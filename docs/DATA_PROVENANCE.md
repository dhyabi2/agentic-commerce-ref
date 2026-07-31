# Data provenance

All downloads performed 2026-07-31 (session date) from NCBI GEO over HTTPS.

## Series retrieved

| Accession | Organism | Samples | Supplementary file | Bytes | What it is |
|---|---|---|---|---|---|
| GSE288028 | Homo sapiens; Mus musculus | 14 | `GSE288028_RAW.tar` (TAR of CellRanger H5) | 302,766,080 | Human pubertal growth plate + 24 h organ culture (vehicle / GH) + 2 mouse growth plate libraries. Chu / Dregval / Trompet / Sävendahl / Chagin, Univ. Gothenburg. Public 2025-12-09. |
| GSE288529 | Mus musculus | 1 | `GSE288529_RAW.tar` (TAR of H5 + MTX/TSV) | 234,659,840 | Single female 4-week-old C57BL/6J mouse, epiphyses digestion. Otsuru lab. PMID 40025030. Public 2025-02-06. |

### sha256

```
3d252c3cd3923651762a31e04479a4ffdbfad4b031c1cb093459c558e7c68c2c  GSE288028_RAW.tar
91aeed9ef961fc505082489c31b0c54844b70892aa6b95c68b724924c36d4439  GSE288529_RAW.tar
```

Raw sequencing reads for the human samples were **not** deposited in GEO — the submitters
withheld them for participant-privacy reasons. Only CellRanger count matrices are public.
This analysis therefore starts from the submitters' count matrices and inherits their
alignment and cell-calling decisions.

## Per-sample libraries (GSE288028)

| GSM | File | Donor/batch | Arm | Species | CellRanger | Genes x cells |
|---|---|---|---|---|---|---|
| GSM9328218 | GSM9328218_P30453_1001.h5 | D1_P30453 | uncultured | human | 7.1.0 | 36601 x 5295 |
| GSM9328219 | GSM9328219_P30453_1002.h5 | D1_P30453 | vehicle | human | 7.1.0 | 36601 x 5133 |
| GSM9328220 | GSM9328220_P30453_1003.h5 | D1_P30453 | GH | human | 7.1.0 | 36601 x 1005 |
| GSM9328221 | GSM9328221_P31011_1001.h5 | D2_P31011 | uncultured | human | 7.1.0 | 36601 x 12911 |
| GSM9328222 | GSM9328222_P31011_1002.h5 | D2_P31011 | vehicle | human | 7.1.0 | 36601 x 3332 |
| GSM9328223 | GSM9328223_P31011_1003.h5 | D2_P31011 | GH | human | 7.1.0 | 36601 x 3840 |
| GSM9328224 | GSM9328224_P25452_001.h5 | D3_P25452 | uncultured | human | 6.0.1 | 36601 x 9115 |
| GSM9328225 | GSM9328225_P25452_004.h5 | D3_P25452 | vehicle (file1) | human | 6.0.1 | 36601 x 10968 |
| GSM9328226 | GSM9328226_P25452_005.h5 | D3_P25452 | vehicle (file2) | human | 6.0.1 | 36601 x 5611 |
| GSM9328227 | GSM9328227_P25452_007.h5 | D3_P25452 | GH (file1) | human | 6.0.1 | 36601 x 7604 |
| GSM9328228 | GSM9328228_P25452_008.h5 | D3_P25452 | GH (file2) | human | 6.0.1 | 36601 x 2868 |
| GSM9328229 | GSM9328229_P22202_1015.h5 | D4_P22202 | uncultured | human | 6.0.1 | 36601 x 383 |
| GSM9328230 | GSM9328230_P27153_1001.h5 | M1_P27153 | uncultured | mouse | 7.0.1 | 32285 x 4712 |
| GSM9328231 | GSM9328231_P27153_1002.h5 | M2_P27153 | uncultured | mouse | 7.0.1 | 32285 x 7157 |

## Age metadata — the decisive provenance fact

GEO records **five** characteristic tags for GSE288028 and no others:
`tissue`, `cell type`, `genotype`, `treatment`, `batch`.

Verified against both the SOFT sample records
(`acc.cgi?acc=GSE288028&targ=gsm&form=text&view=brief`) and the MINiML family XML
(`GSE288028_family.xml`). There is **no** `age`, `sex`, `Tanner stage`, or `development
stage` tag on any sample.

The only age information in the record is a single sentence in the series-level overall
design: *"four human adolescent growth plate chondrocytes were isolated from patients with
age 11-14 years of age undergoing Epiphysiodesis."* That is a group range. It is not
resolved to individual donors, so it cannot be joined to any sample.

(The task brief describes the donors as "aged 12-15, Tanner B2-B4". The deposited GEO
record says 11-14 years and records no Tanner stage. The discrepancy is noted; neither
version supplies per-donor values.)

Donor is additionally confounded with processing batch **and** with CellRanger version
(6.0.1 for D3/D4, 7.1.0 for D1/D2), so even if per-donor ages were recovered from the
authors, an age effect could not be separated from a pipeline-version effect in this design.

GSE288529 is a single animal at a single age (4 weeks), so it contributes no age contrast
either.

## Supplementary rat dataset (analysis E)

| Accession | Organism | File | Bytes | What it is |
|---|---|---|---|---|
| GSE16981 | Rattus norvegicus | `GSE16981_series_matrix.txt.gz` | 2,481,399 | Proximal tibial growth plate, zones separated by manual microdissection, Affymetrix Rat 230 2.0 (GPL1355), MAS5 signal, 35 arrays, n=5 per group, ages 1/3/6/9/12 weeks. Lui, Andrade, Hegde, Chen, Baron, Nilsson. PMIDs 20096814, 20974641. |
| GPL1355 | — | `GPL1355.annot.gz` | 5,714,896 | Platform annotation used for probe -> gene symbol mapping. |

Per-sample characteristics recorded in this series (unlike GSE288028): `age`, `castrated`,
`growth plate zone`. Animals older than 1 week were castrated at day 18.

### sha256 (all raw inputs)

```
3d252c3cd3923651762a31e04479a4ffdbfad4b031c1cb093459c558e7c68c2c  GSE288028_RAW.tar
91aeed9ef961fc505082489c31b0c54844b70892aa6b95c68b724924c36d4439  GSE288529_RAW.tar
959b0339252170160d81d1e4a218c706910d0deef9caf0ef583e57291c316d08  GPL1355.annot.gz
66b8299d510b3fd550793940702d51738ece17b747e7fbeb1f82f5d3e6602eb6  GSE16981_series_matrix.txt.gz
```
