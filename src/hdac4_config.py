"""Fixed configuration for the HDAC4 growth plate analysis.

Everything in this module is pre-registered in docs/PREREGISTRATION.md and must not be
changed after results are inspected.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "results" / "figures"
DOCS = ROOT / "docs"

for _p in (PROC, TABLES, FIGS):
    _p.mkdir(parents=True, exist_ok=True)

SEED = 0

# --------------------------------------------------------------------------------------
# Zone-assignment markers (must not overlap the primary scoring set)
# --------------------------------------------------------------------------------------
ZONE_MARKERS: dict[str, list[str]] = {
    "resting": ["SFRP5", "APOE", "CLU", "PTCH1", "CYTL1", "GAS1", "PTHLH", "RAMP3"],
    "proliferative": ["CCND1", "MKI67", "E2F1"],
    "prehypertrophic": ["IHH", "ALPL"],
    "hypertrophic": ["COL10A1", "MMP13"],
}
ZONE_ORDER = ["resting", "proliferative", "prehypertrophic", "hypertrophic"]

# --------------------------------------------------------------------------------------
# HDAC4 repression-target sets
# --------------------------------------------------------------------------------------
# Full set as given in the brief.
TARGETS_FULL = ["COL10A1", "MMP13", "IBSP", "SPP1", "VEGFA", "ALPL", "PANX3", "SP7"]

# Genes appearing in BOTH the zone-marker list and the target set. Excluded from the
# primary score so that the zonal sanity check is not circular (see PREREGISTRATION 3.1).
TARGETS_OVERLAP = sorted(
    set(TARGETS_FULL) & {g for gs in ZONE_MARKERS.values() for g in gs}
)  # -> ALPL, COL10A1, MMP13

TARGETS_PRIMARY = [g for g in TARGETS_FULL if g not in TARGETS_OVERLAP]
TARGETS_SECONDARY = list(TARGETS_FULL)

# Reported directly, never folded into the score.
PATHWAY_COMPONENTS = [
    "HDAC4", "HDAC5", "HDAC7",
    "SIK1", "SIK2", "SIK3",
    "CAMK4", "CAMK2D",
    "PPP2CA", "PPP2R1A",
    "PTH1R", "PTHLH",
    "MEF2C", "RUNX2",
    "SMAD2", "SMAD3",
    "HIF1A", "EPAS1",
]

# Negative control: identical scoring procedure, must show no zonal gradient.
HOUSEKEEPING = [
    "ACTB", "GAPDH", "TBP", "RPL13A", "PPIA", "SDHA",
    "B2M", "HPRT1", "PGK1", "YWHAZ", "RPLP0", "UBC",
]

# Stress / dissociation signature regressed out during QC.
STRESS_GENES = [
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A", "HSPA1B",
    "HSPB1", "HSP90AA1", "DNAJB1", "SOCS3", "IER2", "ZFP36", "DUSP1", "NR4A1",
]

# Regev-lab cell cycle genes.
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# --------------------------------------------------------------------------------------
# Human -> mouse ortholog symbol mapping (title-case rule plus explicit exceptions)
# --------------------------------------------------------------------------------------
MOUSE_EXCEPTIONS = {"HPRT1": "Hprt", "PIMREG": "Fam64a", "JPT1": "Hn1", "CENPU": "Mlf1ip"}


def to_mouse(sym: str) -> str:
    """Human symbol -> mouse ortholog symbol."""
    if sym in MOUSE_EXCEPTIONS:
        return MOUSE_EXCEPTIONS[sym]
    return sym.capitalize()


def map_genes(symbols: list[str], species: str) -> list[str]:
    return [to_mouse(s) for s in symbols] if species == "mouse" else list(symbols)


# --------------------------------------------------------------------------------------
# GSE288028 sample table.
#
# NOTE: GEO records NO per-donor age, sex or Tanner stage for this series. The only age
# information published is a series-level range ("age 11-14 years") that is not resolved
# to individual donors. `age_years` is therefore None for every donor -- deliberately.
# Donor index is NOT a proxy for age and must never be used as one.
# --------------------------------------------------------------------------------------
GSE288028_SAMPLES = [
    # gsm,          file,                          donor,     arm,        species, cellranger
    ("GSM9328218", "GSM9328218_P30453_1001.h5", "D1_P30453", "uncultured", "human", "7.1.0"),
    ("GSM9328219", "GSM9328219_P30453_1002.h5", "D1_P30453", "vehicle",    "human", "7.1.0"),
    ("GSM9328220", "GSM9328220_P30453_1003.h5", "D1_P30453", "gh",         "human", "7.1.0"),
    ("GSM9328221", "GSM9328221_P31011_1001.h5", "D2_P31011", "uncultured", "human", "7.1.0"),
    ("GSM9328222", "GSM9328222_P31011_1002.h5", "D2_P31011", "vehicle",    "human", "7.1.0"),
    ("GSM9328223", "GSM9328223_P31011_1003.h5", "D2_P31011", "gh",         "human", "7.1.0"),
    ("GSM9328224", "GSM9328224_P25452_001.h5",  "D3_P25452", "uncultured", "human", "6.0.1"),
    ("GSM9328225", "GSM9328225_P25452_004.h5",  "D3_P25452", "vehicle",    "human", "6.0.1"),
    ("GSM9328226", "GSM9328226_P25452_005.h5",  "D3_P25452", "vehicle",    "human", "6.0.1"),
    ("GSM9328227", "GSM9328227_P25452_007.h5",  "D3_P25452", "gh",         "human", "6.0.1"),
    ("GSM9328228", "GSM9328228_P25452_008.h5",  "D3_P25452", "gh",         "human", "6.0.1"),
    ("GSM9328229", "GSM9328229_P22202_1015.h5", "D4_P22202", "uncultured", "human", "6.0.1"),
    ("GSM9328230", "GSM9328230_P27153_1001.h5", "M1_P27153", "uncultured", "mouse", "7.0.1"),
    ("GSM9328231", "GSM9328231_P27153_1002.h5", "M2_P27153", "uncultured", "mouse", "7.0.1"),
]

# Per-donor age: unavailable. Recorded explicitly so downstream code fails loudly rather
# than silently substituting a surrogate.
DONOR_AGE_YEARS: dict[str, float | None] = {
    "D1_P30453": None,
    "D2_P31011": None,
    "D3_P25452": None,
    "D4_P22202": None,
    "M1_P27153": None,
    "M2_P27153": None,
}
AGE_RANGE_SERIES_LEVEL = "11-14 years (four donors, not resolved per donor)"

# QC thresholds
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3
MAX_PCT_MT = 20.0
MIN_COUNTS_PER_CELL = 500
MIN_CELLS_PER_ZONE = 100  # pre-registered flag threshold
N_PERMUTATIONS = 1000
