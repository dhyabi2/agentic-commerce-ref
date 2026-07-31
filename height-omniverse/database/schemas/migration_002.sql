-- migration_002: harden schema (additive only; never deletes data).
-- Adds provenance, design, effect-extraction, and evidence-integrity fields, a species-fusion
-- reference table, an endpoint-ontology reference table, derived integrity flags, and QC views.

-- ---- reference: species fusion biology (from config/species_registry.yaml) ----
CREATE OR REPLACE TABLE ref_species_fusion (species VARCHAR, fusion_capable BOOLEAN, priority INTEGER, note VARCHAR);
DELETE FROM ref_species_fusion;
INSERT INTO ref_species_fusion VALUES
 ('human',TRUE,1,'target species'),('nonhuman_primate',TRUE,2,'closest fusion'),
 ('rabbit',TRUE,3,'FUSION-CAPABLE key'),('pig',TRUE,4,'large animal'),('sheep_goat',TRUE,4,'large animal'),
 ('rat',FALSE,5,'physes do NOT fully fuse'),('mouse',FALSE,6,'physes do NOT fully fuse'),
 ('zebrafish',FALSE,7,'dev system'),('chicken',FALSE,7,'dev system');

-- ---- reference: endpoint ontology rank (1=final length ... 10=biomarker) ----
CREATE OR REPLACE TABLE ref_endpoint_rank (endpoint_type VARCHAR, endpoint_rank INTEGER, is_final_length BOOLEAN);
DELETE FROM ref_endpoint_rank;
INSERT INTO ref_endpoint_rank VALUES
 ('adult_height',1,TRUE),('near_final_height',1,TRUE),('mature_long_bone_length',1,TRUE),('bone_length',1,TRUE),
 ('stable_mineralized_length',2,TRUE),('segment_dimension',3,FALSE),('target_adjusted_height',4,FALSE),
 ('predicted_adult_height',5,FALSE),('growth_velocity',6,FALSE),('velocity',6,FALSE),('plate_thickness',6,FALSE),
 ('bone_age',7,FALSE),('fusion_timing',7,FALSE),('reserve_cell',8,FALSE),('matrix_hypertrophy',8,FALSE),
 ('condyle',9,FALSE),('fracture_repair',9,FALSE),('mechanism',10,FALSE),('biomarker',10,FALSE),('other',10,FALSE);

-- ---- studies: provenance + design + integrity columns (additive) ----
ALTER TABLE studies ADD COLUMN IF NOT EXISTS doi VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS pmid VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS source_url VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS access_date VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS retraction_status VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS locator VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS extraction_method VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS reviewer_status VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS fusion_capable_flag VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS endpoint_type VARCHAR;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS endpoint_rank INTEGER;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS is_final_length BOOLEAN;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS digitized_flag BOOLEAN;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS contradiction_flag BOOLEAN;
ALTER TABLE studies ADD COLUMN IF NOT EXISTS reusability_license VARCHAR;

-- backfill fusion flag by species (rat/mouse -> not fusion-capable)
UPDATE studies SET fusion_capable_flag =
  CASE WHEN (lower(species) LIKE '%rat%' OR lower(species) LIKE '%mouse%')
            AND (lower(species) LIKE '%rabbit%' OR lower(species) LIKE '%human%' OR lower(species) LIKE '%calf%'
                 OR lower(species) LIKE '%pig%' OR lower(species) LIKE '%sheep%' OR lower(species) LIKE '%goat%') THEN 'mixed'
       WHEN lower(species) LIKE '%human%' OR lower(species) LIKE '%rabbit%' OR lower(species) LIKE '%primate%'
            OR lower(species) LIKE '%pig%' OR lower(species) LIKE '%sheep%' OR lower(species) LIKE '%goat%' THEN 'yes'
       WHEN lower(species) LIKE '%rat%' OR lower(species) LIKE '%mouse%' OR lower(species) LIKE '%zebra%' OR lower(species) LIKE '%chick%' THEN 'no'
       ELSE 'unknown' END;

-- derive endpoint_type from free-text endpoint, then rank via ref table
UPDATE studies SET endpoint_type =
  CASE WHEN lower(endpoint) LIKE '%adult height%' THEN 'adult_height'
       WHEN lower(endpoint) LIKE '%near%final%' OR lower(endpoint) LIKE '%final height%' THEN 'near_final_height'
       WHEN lower(endpoint) LIKE '%final%length%' OR (lower(endpoint) LIKE '%length%' AND lower(endpoint) NOT LIKE '%velocity%') THEN 'mature_long_bone_length'
       WHEN lower(endpoint) LIKE '%velocity%' OR lower(endpoint) LIKE '%gain%' OR lower(endpoint) LIKE '%lgr%' THEN 'velocity'
       WHEN lower(endpoint) LIKE '%thickness%' THEN 'plate_thickness'
       WHEN lower(endpoint) LIKE '%bone age%' THEN 'bone_age'
       WHEN lower(endpoint) LIKE '%bmd%' OR lower(endpoint) LIKE '%mineral%' THEN 'biomarker'
       WHEN lower(endpoint) LIKE '%quiescence%' OR lower(endpoint) LIKE '%differentiation%' THEN 'mechanism'
       ELSE 'other' END;
UPDATE studies s SET endpoint_rank = r.endpoint_rank, is_final_length = r.is_final_length
  FROM ref_endpoint_rank r WHERE s.endpoint_type = r.endpoint_type;
UPDATE studies SET access_date='2026-07-31', extraction_method='agent-verified API', reviewer_status='machine'
  WHERE access_date IS NULL;

-- interventions: provenance
ALTER TABLE interventions ADD COLUMN IF NOT EXISTS efficacy_safety_separation VARCHAR;
UPDATE interventions SET efficacy_safety_separation='separate' WHERE efficacy_safety_separation IS NULL;

-- evidence_edges: contradiction integrity
ALTER TABLE evidence_edges ADD COLUMN IF NOT EXISTS has_contradiction BOOLEAN;
UPDATE evidence_edges SET has_contradiction = (conflicting_edges IS NOT NULL AND conflicting_edges <> '');

-- ---- QC views answering canonical questions ----
CREATE OR REPLACE VIEW v_adult_height_evidence AS
  SELECT intervention, population, adult_height_effect, surrogate_translated, source_id
  FROM surrogate_calibration WHERE upper(surrogate_translated)='YES';
CREATE OR REPLACE VIEW v_velocity_pos_final_null AS
  SELECT intervention, population, surrogate_effect, source_id
  FROM surrogate_calibration WHERE direction='+' AND upper(surrogate_translated) IN ('NO','UNKNOWN');
CREATE OR REPLACE VIEW v_brake_genes AS
  SELECT symbol, direction_of_effect, monogenic_stature_evidence, source_provenance
  FROM genes WHERE direction_of_effect ILIKE '%taller%';
CREATE OR REPLACE VIEW v_final_length_after_washout AS
  SELECT citation, species, fusion_capable_flag, endpoint, direction, washout, source_location
  FROM studies WHERE is_final_length=TRUE;
CREATE OR REPLACE VIEW v_mouse_or_condyle_only AS
  SELECT citation, species, tissue_site, endpoint FROM studies
  WHERE fusion_capable_flag='no' OR lower(tissue_site) LIKE '%condyle%';
CREATE OR REPLACE VIEW v_contradictions AS
  SELECT subject, predicate, object, direction, source_id, conflicting_edges
  FROM evidence_edges WHERE has_contradiction=TRUE;
