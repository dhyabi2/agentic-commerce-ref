#!/usr/bin/env python
"""Emit curated verified-evidence CSVs into data/curated/ (loaded by build_db.py).
Every row traces to a source_id. Values behind paywall = 'NOT RETRIEVED'. Data verified by
harvest agents via NCBI E-utilities / Europe PMC / ClinicalTrials.gov / GEO (retrieved 2026-07-31)."""
import csv, pathlib
CUR = pathlib.Path(__file__).resolve().parents[2] / "data" / "curated"; CUR.mkdir(parents=True, exist_ok=True)
def write(name, cols, rows):
    with open(CUR / name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(cols); w.writerows(rows)
    print(f"{name}: {len(rows)} rows")

# ---------------- surrogate -> final-height calibration (crown jewel) ----------------
CAL_COLS = ["intervention","population","design","n","surrogate_effect","adult_height_effect","surrogate_translated","direction","source_id","notes"]
CAL = [
["rGH","idiopathic short stature","double-blind RCT vs placebo","68 enrolled/33 to AH","velocity/PAH not separately tabulated","+0.51 SDS (3.7 cm) vs placebo","YES","+","PMID 15240584","Leschek 2004; modest but real translation"],
["rGH","GH deficiency","registry/cohort to final height (KIGS)","369","large velocity gain","final H -1.5 SDS; gain +1.7 SDS; ~87% normalize","YES","+","PMID 10102057","Cutfield/KIGS 1999; reaches target range"],
["Long-acting GH","pediatric GHD","RCT/NMA 52-wk","NA","AHV/HtSDS ~ daily GH","NOT RETRIEVED (no final-height yet)","UNKNOWN","0","PMC11480202","approved; final-height not yet published"],
["Letrozole (aromatase inhibitor)","short/ISS boys","double-blind RCT 2 yr","31","PAH +5.9 cm; HtSDS-for-BA +0.7","NOT MEASURED (PAH only)","UNKNOWN","+","PMID 16189252","Hero 2005; PAH surrogate only"],
["Letrozole (aromatase inhibitor)","ISS boys","RCT to ADULT height","20","(PAH gain at baseline)","AH 164.8 vs 163.7 cm; P=0.49 NS","NO","0","PMID 31024444","Varimo/Hero 2019; PAH gain did NOT translate; vertebral deformities noted"],
["GnRH agonist","short adolescents, normal puberty (ISS-like)","RCT to adult height","50","delays BA/puberty","+4.2 cm/+0.6 SD vs placebo","YES","+","PMID 12621135","Yanovski 2003; BUT decreased BMD; not routinely recommended"],
["GnRH agonist","early-onset central precocious puberty","cohorts/observational","varies","halts BA advance; PAH rises","AH > initial PAH by ~3.9 cm; reaches target","YES","+","PMC5752503","benefit mainly EARLY-onset CPP; little/none if late"],
["GnRH agonist + GH","early-puberty girls, short PAH","RCTs/cohorts","varies","additive PAH/velocity","modest AH gain over GnRHa alone; inconsistent","PARTIAL","+","PMC6332953","combination benefit context-dependent; not established in normal-puberty ISS"],
["Oxandrolone (+GH)","Turner syndrome","double-blind RCT to final height","~106; meta 498","velocity/BA increase","+~4.6 cm final height vs placebo (on GH)","YES","+","PMID 21493672","Gault 2011; late pubertal induction added ~3.8 cm"],
["Vosoritide","achondroplasia","phase-3 double-blind RCT 52 wk","121","AGV +1.57 cm/yr (P<0.0001)","NOT in RCT; ~13.6 cm/8 yr vs natural history (press, non-peer-reviewed)","UNKNOWN","+","PMID 32891212","approved; no controlled adult-height endpoint published"],
["Navepegritide (TransCon CNP)","achondroplasia","phase-3 RCT (ApproaCH) 52 wk","84","AGV 5.84 vs 3.88 cm/yr","NOT RETRIEVED (none)","UNKNOWN","+","PMID 41247754","investigational; adult height not assessed"],
["Infigratinib (FGFR3i)","achondroplasia","Ph2 dose-finding + Ph3 RCT","Ph3 ~110","Ph3 AHV +2.10 cm/yr vs placebo","NOT RETRIEVED (none)","UNKNOWN","+","PMID 42370681","PROPEL; velocity only; highest RCT AHV delta to date"],
["TYRA-300 (FGFR3i)","achondroplasia","Ph2 open-label (BEACH301)","dosing from 2025-08","NOT RETRIEVED (no data yet)","NOT RETRIEVED","UNKNOWN","0","NCT BEACH301","first child dosed 2025-08-21; readout expected 2H2026"],
["Mecasermin (rhIGF-1)","severe primary IGF-1 deficiency / Laron","registry to near-adult height","102","velocity gain","NAH HtSDS gain +0.9 (prepubertal-naive +1.4); Laron 10.5% normal","PARTIAL","+","PMID 40626687","disease context; residual deficit esp. Laron"],
["Estrogen (low-dose childhood, +/-GH)","Turner syndrome","double-blind 2x2 RCT to adult height","~149","modest velocity","GH +0.78 SDS (5.0 cm); estrogen-alone small; GH+estrogen best","YES","+","PMID 21449786","Ross 2011 NEJM; ultra-low-dose estrogen adds alongside GH"],
["Levothyroxine","congenital hypothyroidism (screened)","observational cohorts","varies","catch-up velocity","normal adult height if treated early; incomplete if late","YES","+","PMID 11580992","target/genetic height dominant; acquired hypoT catch-up may be incomplete"],
["Nutrition/energy supplementation","stunting/undernutrition (LMIC)","RCTs/systematic reviews","many","short-term catch-up, attenuates after stopping","NO robust final-height recovery after ~2 yr","NO","+","PMC7360439","window-limited; prevention > reversal"],
]
write("surrogate_calibration.csv", CAL_COLS, CAL)

# ---------------- interventions ----------------
IV_COLS = ["intervention_id","name","synonyms","modality","target","action","approval_status","route","pk_summary","evidence_context","technology_readiness","source_provenance"]
IV = [
["IV_GH","Recombinant GH (somatropin)","somatropin","hormone (recombinant protein)","GH receptor","agonist (JAK/STAT->IGF-1)","approved (ISS/GHD/Turner)","SC daily","short t1/2 (hrs; NOT RETRIEVED)","human final-height RCT/registry","approved","PMID 15240584; 10102057"],
["IV_LAGH","Long-acting GH","lonapegsomatropin; somatrogon; somapacitan","hormone (modified GH)","GH receptor","agonist","approved (pediatric GHD)","SC weekly","weekly; t1/2 NOT RETRIEVED","52-wk velocity RCTs","approved","PMC11480202"],
["IV_AI","Aromatase inhibitor","letrozole; anastrozole","small molecule","aromatase (CYP19)","inhibitor (blocks estrogen synthesis)","approved oncology; OFF-LABEL short stature","oral daily","letrozole t1/2 ~2 d (UNVERIFIED)","PAH RCT + adult-height RCT","approved (off-label use)","PMID 16189252; 31024444"],
["IV_GNRHA","GnRH agonist","leuprolide; triptorelin; histrelin","peptide analogue","GnRH receptor","agonist (pituitary desensitization)","approved (CPP)","SC/IM depot or implant","depot 1-6 mo","adult-height RCT (CPP/ISS)","approved","PMID 12621135; PMC5752503"],
["IV_OXA","Oxandrolone","","small molecule (anabolic androgen)","androgen receptor","agonist","approved (adjunct; Turner)","oral daily","t1/2 NOT RETRIEVED","final-height RCT (Turner)","approved","PMID 21493672"],
["IV_VOSO","Vosoritide","BMN-111","peptide (CNP analogue)","NPR-B (NPR2)","agonist (inhibits FGFR3/MAPK)","approved (achondroplasia)","SC daily","short t1/2 (engineered)","phase-3 velocity RCT","approved","PMID 32891212"],
["IV_NAVE","Navepegritide","TransCon CNP","peptide prodrug (PEG sustained release)","NPR-B (NPR2)","agonist","investigational (BLA-stage)","SC weekly","weekly release","phase-3 velocity RCT","phase 3","PMID 41247754"],
["IV_INFI","Infigratinib","BGJ398; Truseltiq","small-molecule TKI","FGFR1-3","inhibitor","investigational (onc-approved)","oral daily","t1/2 NOT RETRIEVED","phase-3 velocity RCT","phase 3","PMID 42370681"],
["IV_TYRA","TYRA-300","dabogratinib","small-molecule TKI","FGFR3-selective","inhibitor","investigational (Ph2)","oral","NOT RETRIEVED","no efficacy data yet","phase 2","NCT BEACH301"],
["IV_IGF1","Mecasermin (rhIGF-1)","","peptide hormone","IGF-1 receptor","agonist","approved (severe primary IGF-1 deficiency)","SC BID","short t1/2 (BID)","near-adult-height registry","approved","PMID 40626687"],
["IV_ESTR","Estrogen (low-dose)","ethinylestradiol; estradiol","steroid hormone","estrogen receptor","agonist (growth-plate + GH/IGF axis)","approved (Turner induction)","oral/transdermal","NOT RETRIEVED","adult-height RCT (Turner)","approved","PMID 21449786"],
["IV_LT4","Levothyroxine","","thyroid hormone (T4)","thyroid hormone receptor","agonist","approved","oral daily","t1/2 ~7 d (UNVERIFIED)","cohort to adult height","approved","PMID 11580992"],
["IV_NUTR","Nutrition/energy supplementation","","nutrient/public-health","systemic anabolic/IGF-1 substrate","substrate","NA","oral","NA","RCT/reviews (stunting)","established","PMC7360439"],
["IV_ABL","Abaloparatide","BA058; TYMLOS","peptide (PTHrP analogue)","PTH1R","agonist (RG-biased, transient cAMP)","approved (osteoporosis)","SC daily","Tmax 0.51 h; t1/2 ~1 h","condyle only in growing organisms; NO conventional-physis final length","approved (diff. indication)","companion analysis; PMID 26562265"],
]
write("interventions.csv", IV_COLS, IV)

# ---------------- evidence edges: intervention -> attained adult height ----------------
E_COLS = ["edge_id","subject","predicate","object","direction","context","effect_size","source_id","evidence_tier","directness","confidence","conflicting_edges"]
E = [
["E1","Recombinant GH","changes","attained adult height","+","ISS, RCT","+0.51 SDS (3.7 cm)","PMID 15240584","T1","high","moderate",""],
["E2","Recombinant GH","changes","attained adult height","+","GHD, registry","+1.7 SDS","PMID 10102057","T1","high","high",""],
["E3","Aromatase inhibitor (letrozole)","changes","predicted adult height","+","ISS boys, RCT","PAH +5.9 cm","PMID 16189252","T1","high","moderate","E4"],
["E4","Aromatase inhibitor (letrozole)","changes","attained adult height","0","ISS boys, RCT to AH","NS (P=0.49)","PMID 31024444","T1","high","moderate","E3"],
["E5","GnRH agonist","changes","attained adult height","+","early-onset CPP","AH > initial PAH ~3.9 cm","PMC5752503","T1","med","moderate","E6"],
["E6","GnRH agonist","changes","attained adult height","0","normal-puberty ISS","small; not recommended","PMID 12621135","T1","high","moderate","E5"],
["E7","Oxandrolone (+GH)","changes","attained adult height","+","Turner, RCT","+4.6 cm","PMID 21493672","T1","high","moderate",""],
["E8","Vosoritide","changes","annualized growth velocity","+","achondroplasia, RCT","+1.57 cm/yr","PMID 32891212","T3","high","high","E9"],
["E9","Vosoritide","changes","attained adult height","?","achondroplasia","no controlled endpoint","PMID 32891212","T8","low","low","E8"],
["E10","Infigratinib","changes","annualized growth velocity","+","achondroplasia, RCT","+2.10 cm/yr","PMID 42370681","T3","high","moderate",""],
["E11","Navepegritide","changes","annualized growth velocity","+","achondroplasia, RCT","+1.96 cm/yr","PMID 41247754","T3","high","moderate",""],
["E12","Mecasermin (rhIGF-1)","changes","near-final height","+","primary IGF-1 deficiency","+0.9 SDS (partial)","PMID 40626687","T1","med","moderate",""],
["E13","Estrogen (low-dose)","changes","attained adult height","+","Turner, RCT (adjunct)","adjunct to GH","PMID 21449786","T1","high","moderate",""],
["E14","Levothyroxine","changes","attained adult height","+","early-treated hypothyroidism","normal if early","PMID 11580992","T1","high","moderate",""],
["E15","Nutrition/energy supplementation","changes","attained adult height","0","stunting after ~2 yr","no robust recovery","PMC7360439","T1","med","moderate",""],
["E16","Abaloparatide","changes","mandibular condyle length","+","adolescent rat, condyle (secondary cartilage)","+1.34 mm (local)","PMID 37727017","T6","low","moderate",""],
["E17","Abaloparatide","changes","conventional long-bone final length","?","no study exists","not identified","companion analysis","T8","none","low",""],
]
write("evidence_edges.csv", E_COLS, E)

# ---------------- access failures (blocked/unretrieved) ----------------
AF_COLS = ["source","reason","attempted_method","metadata_found","request_route","importance"]
AF = [
["Long-acting GH final-height data","not yet published","Europe PMC/registry","52-wk velocity only","await publication","high"],
["Vosoritide controlled adult-height endpoint","does not exist (velocity-only trials)","PubMed/registry","press ~13.6 cm/8 yr vs nat-hist","await BioMarin LTE peer review","high"],
["Navepegritide/Infigratinib/TYRA-300 adult height","velocity-only / no data yet","PubMed/ClinicalTrials","AGV deltas","await final-height endpoints","high"],
["Letrozole 2020-2025 final-height RCT (>N=20)","not located","Europe PMC","only Varimo 2019 N=20","author/registry","medium"],
["Precise PK half-lives (most compounds)","FDA-label PK pages not fetched","label lookup","dosing frequency only","fetch DailyMed/label","medium"],
["Restricted individual-level genetics (UKB/dbGaP/EGA)","access-controlled","N/A (public summary only)","summary stats only","formal data access application","medium"],
]
AF += [
["PEMF anchor (Yang/Sung 2021)","not MEDLINE-indexed (Applied Sciences, MDPI)","Crossref/DOI","GP length +25.6%, IGF-1 +13.6%; DOI 10.3390/app11167571","cite DOI; endpoint is GP thickness NOT final length","medium"],
["Final adult bone-length effect of ANY mechanical intervention in humans","does not exist","Europe PMC","only velocity/thickness or no-effect (LIPUS null)","decisive experiment needed","high"],
["Botanical final-height evidence (icariin, Cissus, Eucommia, Phlomis, p-coumaric)","only animal velocity / none for longitudinal growth","Europe PMC","HT042 human RCT +0.29 cm 24-wk velocity only (PMID 29130588)","human final-height RCTs","medium"],
["Yengo 2022 growth-plate/chondrocyte cell-type enrichment","not in retrieved abstract","Europe PMC","abstract says only 'biologically relevant genes'","full-text/supplement + Wood 2014","medium"],
["GSE288028 human raw sequence","withheld (PII)","GEO","processed matrices public only","controlled-access request","medium"],
]
write("access_failures.csv", AF_COLS, AF)

# ---------------- genetic direction of effect (verified PMIDs) ----------------
GD_COLS = ["gene","variant_class","activity_change","height_direction","dev_vs_postnatal","pathway","key_phenotype","source_id"]
GD = [
["SHOX","LoF/haploinsufficiency","down","shorter","developmental+postnatal","growth-plate homeobox TF","Leri-Weill dyschondrosteosis / ISS; mesomelic","PMID 20301394"],
["ACAN","LoF (het)","down","shorter","growth-plate; postnatal cessation","aggrecan matrix","short + ADVANCED bone age, early OA","PMID 27870580"],
["NPR2","LoF (biallelic)","down","shorter","developmental+postnatal","CNP->NPR2/cGMP","acromesomelic dysplasia Maroteaux; het short","PMID 15146390"],
["NPR2","GoF","up","taller","developmental+postnatal","excess cGMP","skeletal overgrowth/tall","PMID 22870295"],
["NPPC","overexpression/dup","up","taller","developmental+postnatal","CNP ligand -> NPR2","overgrowth, Marfanoid","PMID 17373680"],
["NPR3","LoF (biallelic)","up (reduced CNP clearance)","taller","developmental+postnatal","loss of CNP clearance receptor","tall + scoliosis","PMID 40171685"],
["FGFR3","GoF (G380R)","up (constitutive)","shorter","developmental+postnatal","FGFR3 inhibits proliferation","achondroplasia","PMID 7913883"],
["FGFR3","GoF (N540K)","up","shorter","developmental+postnatal","FGFR3","hypochondroplasia","PMID 8589686"],
["FGFR3","LoF (R621H)","down","taller","developmental+postnatal","loss of growth-plate braking","CATSHL: tall stature","PMID 17033969"],
["IHH","LoF (missense)","down","shorter (digits)","developmental (digit patterning)","Indian hedgehog","brachydactyly type A1","PMID 20024692"],
["GHR","LoF","down","shorter","postnatal","GH receptor (systemic)","Laron syndrome; low IGF-1","PMID 8488849"],
["STAT5B","LoF","down","shorter","postnatal","GHR->STAT5B->IGF1","GH insensitivity + immunodeficiency","PMID 13679528"],
["IGF1","LoF (deletion)","down","shorter","prenatal+postnatal","IGF-1 ligand (systemic)","IUGR + postnatal failure, deafness","PMID 8857020"],
["IGF1R","LoF (het)","down","shorter","prenatal+postnatal","IGF-1 receptor (systemic)","SGA/IUGR, postnatal short","PMID 14657428"],
["IGFALS","LoF","down","shorter (mild)","postnatal","acid-labile subunit","ALS deficiency; very low IGF-1/IGFBP-3","PMID 21396577"],
["PAPPA2","LoF","down (reduced free IGF-1)","shorter","postnatal","protease liberating IGF-1 from IGFBP","progressive growth failure; high total/low free IGF-1","PMID 26902202"],
["STC2","rare coding (inhibits PAPP-A)","up (more free IGF-1)","taller","postnatal","stanniocalcin-2 inhibits PAPP-A","height-increasing coding variant","PMID 28146470"],
["CYP19A1","LoF (aromatase deficiency)","down (estrogen synthesis)","taller","postnatal (fusion boundary)","estrogen biosynthesis (systemic)","tall, UNFUSED epiphyses, continued growth","PMID 8530621"],
["ESR1","LoF (estrogen resistance)","down (estrogen action)","taller","postnatal (fusion boundary)","ERalpha in growth plate","tall (204 cm), unfused epiphyses","PMID 8090165"],
["DNMT3A","LoF (germline het)","down (methylation)","taller","developmental+postnatal","DNA methyltransferase","Tatton-Brown-Rahman overgrowth","PMID 27701732"],
["NSD1","LoF/haploinsufficiency","down","taller","prenatal+postnatal","histone methyltransferase","Sotos overgrowth","PMID 12687502"],
["EZH2","LoF (het)","down","taller","prenatal+postnatal","PRC2 H3K27 methyltransferase","Weaver overgrowth","PMID 24214728"],
["GPC3","LoF (X-linked)","down","taller","prenatal+postnatal","glypican-3 (IGF2/hedgehog)","Simpson-Golabi-Behmel overgrowth","PMID 8589713"],
["IGF2","LoF (paternal)","down","shorter","prenatal (imprinted)","imprinted IGF-2 ligand","Silver-Russell growth restriction","PMID 26154720"],
]
write("genetic_direction.csv", GD_COLS, GD)

# genetics edges (gene activity -> adult height direction); map taller/shorter -> +/-
def _dir(h):
    h = h.lower()
    return "+" if "taller" in h else "-" if "short" in h else "?"
GE = [[f"G{i+1}", r[0], "altered-activity-changes", "adult height", _dir(r[3]), f"{r[1]}; {r[4]}", r[6], r[7], "T5", "high (human genetics)", "moderate", ""]
      for i, r in enumerate(GD)]
# ---------------- studies: mechanical + nutrient + selected genetics ----------------
S_COLS = ["study_id","citation","year","design","species","sex","age_stage","disease","tissue_site","intervention","comparator","n_true_bioreplicates","exposure","duration","washout","endpoint","maturity_followup","effect_size","direction","geometry_mineralization","directness_tier","risk_of_bias","source_location","notes"]
S = [
# study_id,citation,year,design,species,sex,age_stage,disease,tissue_site,intervention,comparator,n_true,exposure,duration,washout,endpoint,maturity_followup,effect_size,direction,geometry_mineralization,directness_tier,risk_of_bias,source_location,notes
["S_PEMF","Yang/Sung 2021 Applied Sciences (DOI 10.3390/app11167571)",2021,"animal experiment","rat (immature)","NR","3 wk","none","proximal tibia growth plate","PEMF 28Hz 20G/2mT 50% duty 10h/d","sham","NR","28Hz/2mT/50%/10h-d","10 d","N","GP length; serum IGF-1","N","GP length +25.6%; IGF-1 +13.6%","+","N","T3","moderate","DOI 10.3390/app11167571","GP thickness/velocity proxy NOT final length; not MEDLINE-indexed"],
["S_LIPUS","Spadaro & Albanese 1998",1998,"animal experiment","rat (young)","NR","juvenile","none","femur/tibia","LIPUS 1.5MHz 30mW/cm2 20min/d","control limb","NR","20min/day","28 d","N","FINAL femur/tibia length; BMD","N","no effect (<=2% length)","0","Y","T2","low","PMID 9651966","final length measured -> NULL"],
["S_HV","Stokes 2006",2006,"animal experiment","rat+rabbit+calf (immature)","NR","immature","none","growth plate","static compression vs distraction","contralateral","NR","sustained","8 d","N","longitudinal growth rate (fluorochrome)","N","static SLOWS; distraction ACCELERATES","+/-","N","T3","low","PMID 16705695","Hueter-Volkmann; mechanical direction; mixed species"],
["S_JUMP","Jumping exercise RCT 2025",2025,"RCT","human","NR","8-11 y (prepubertal, short)","short stature","whole body","progressive jumping 3x/wk","control","NR","50min 3x/wk","24 wk","N","24-wk height gain; femoral-neck BMD","N","+4.20 vs 2.48 cm (p=0.001)","+","Y","T3","moderate","PMID 41233903","VELOCITY over 24 wk, NOT final height"],
["S_LOAD","McGarry/De Luca 2024",2024,"animal experiment","mouse","F","juvenile","none","tibial growth plate","cyclic physiologic loading","contralateral","NR","cyclic","4 wk","N","longitudinal growth; GP height; PTHrP","N","increased (PTHrP-mediated)","+","N","T3","low","PMID 39090666","endogenous PTHrP mechanism"],
["S_HT042","HT042 pediatric RCT 2018",2018,"RCT double-blind","human","M/F","<25th pctile children","short stature","systemic","HT042 (Astragalus mix)","placebo","129","oral","24 wk","N","24-wk height gain; bone age; IGF-1","N","+0.29 cm vs placebo (p=0.027)","+","N","T3","moderate","PMID 29130588","botanical; velocity not final height; larger effect <10th pctile"],
["S_ZINC","Imdad zinc meta-analysis 2011",2011,"meta-analysis of RCTs","human","M/F","<5 y (LMIC)","deficiency risk","systemic","zinc supplementation","placebo","many","~10 mg/d","varies","N","linear growth (HAZ/length)","N","+0.13 SMD (~+0.37 cm/24wk)","+","N","T1","low","PMID 21501440","deficiency-correction; small"],
["S_VITD","Vitamin D linear-growth meta 2020",2020,"Cochrane meta","human","M/F","<5 y","varies","systemic","vitamin D","placebo","many","varies","varies","N","linear growth","N","little-to-no effect (low certainty)","0","N","T1","low","PMID 33305842","replete populations no benefit"],
["S_PROT","Kristjansson food-supplementation meta 2015",2015,"Cochrane meta","human","M/F","3mo-5y (disadvantaged)","undernutrition","systemic","protein/energy supplementation","control","many","food","varies","N","growth","N","small increase","+","N","T1","low","PMID 25739460","window-limited"],
["S_ACAN","Gkourogianni ACAN cohort 2017",2017,"case series (genetics)","human","M/F","children->adult","ACAN LoF","growth plate","(natural experiment)","unaffected","103","germline","lifelong","NA","adult height; bone age","Y","short + advanced bone age","-","Y","T5","moderate","PMID 27870580","allelic series boundary"],
["S_CYP19","Morishima aromatase deficiency 1995",1995,"case report (genetics)","human","M","adult","CYP19A1 LoF","fusion boundary","(natural experiment)","norms","1","germline","lifelong","NA","adult height; epiphyseal fusion","Y","tall; UNFUSED epiphyses; continued growth","+","Y","T5","moderate","PMID 8530621","estrogen governs fusion"],
]
assert all(len(r)==len(S_COLS) for r in S), [ (r[0],len(r)) for r in S if len(r)!=len(S_COLS)]
write("studies.csv", S_COLS, S)

# append genetics + mechanical/nutrient edges to the edges file
ME = [
["M1","PEMF (28Hz/2mT)","changes","growth-plate thickness","+","immature rat","+25.6%","DOI 10.3390/app11167571","T3","med","low","M2"],
["M2","LIPUS","changes","final long-bone length","0","young rat","no effect","PMID 9651966","T2","high","moderate","M1"],
["M3","Static compression","changes","longitudinal growth rate","-","immature (Hueter-Volkmann)","slows","PMID 16705695","T3","high","moderate",""],
["M4","Distraction/tension","changes","longitudinal growth rate","+","immature","accelerates","PMID 16705695","T3","high","moderate",""],
["M5","Jumping exercise","changes","24-wk height gain (velocity)","+","prepubertal children","+1.7 cm/24wk","PMID 41233903","T3","med","moderate",""],
["M6","HT042 (botanical)","changes","24-wk height gain (velocity)","+","short children RCT","+0.29 cm","PMID 29130588","T3","med","moderate",""],
["M7","Zinc","changes","linear growth","+","deficiency-risk children","+0.13 SMD","PMID 21501440","T1","med","low",""],
]
with open(CUR / "evidence_edges.csv", "a", newline="") as f:
    w = csv.writer(f); w.writerows(GE); w.writerows(ME)
print(f"appended {len(GE)} genetics + {len(ME)} mechanical/nutrient edges")

# ---------------- standalone human-readable result tables ----------------
RES = pathlib.Path(__file__).resolve().parents[2] / "results" / "tables"; RES.mkdir(parents=True, exist_ok=True)
OMICS_COLS = ["accession","organism","assay","processed_downloadable","url","relevance"]
OMICS = [
["GSE288529","mouse","scRNA-seq 10x","YES","geo/acc?acc=GSE288529","4-wk epiphysis; ApoE resting-zone marker (HIGH)"],
["GSE244884/880/881","mouse","scRNA-seq","YES","geo/acc?acc=GSE244884","Hedgehog/Ptch1 resting-zone, P36, PTHrP-creER lineage (HIGH)"],
["GSE275171","mouse","scRNA-seq","YES","geo/acc?acc=GSE275171","Fgfr3 achondroplasia resting-zone/CREB, P19-22 (HIGH)"],
["GSE288028","human+mouse","scRNA-seq","PARTIAL (human raw withheld PII)","geo/acc?acc=GSE288028","2026 human pubertal GP 11-14y + GH stim (HIGH)"],
["PRJNA793434","mouse","LMD RNA-seq (zonal)","SRA raw","bioproject/PRJNA793434","resting vs proliferative laser-microdissection (HIGH)"],
["PRJNA937633 / GSE225796","mouse","zonal RNA-seq","YES","geo/acc?acc=GSE225796","round/resting-layer specificity + height-GWAS enrichment (HIGH)"],
["PRJNA832296/832298","mouse","scRNA/bulk","SRA raw","bioproject/PRJNA832296","Idh1-mutant EMBRYONIC E18.5 (keep separate)"],
["GSE270640/270641","mouse","expr+methylation","YES","geo/acc?acc=GSE270640","Dnmt1 determines bone length via GP metabolism (HIGH)"],
["GSE160364","mouse","bulk RNA-seq","YES","geo/acc?acc=GSE160364","label-retaining (resting) chondrocytes"],
["GSE271634 / GSE326052 / GSE275252 / GSE182540","mouse","scRNA-seq","YES","GEO","additional postnatal GP scRNA (aging, GPR126, chondro-HSC)"],
["GSE267139","human","scRNA-seq","YES","geo/acc?acc=GSE267139","human cartilage development (largely fetal)"],
]
with open(RES / "OMICS_AVAILABILITY.csv","w",newline="") as f:
    csv.writer(f).writerows([OMICS_COLS]+OMICS)
print(f"OMICS_AVAILABILITY.csv: {len(OMICS)} datasets")
print("curated CSVs written to", CUR)
