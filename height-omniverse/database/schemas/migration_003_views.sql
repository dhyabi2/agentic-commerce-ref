-- migration_003: additional canonical query views (Phase M) + a unified text-search surface.
-- Every view resolves to source-linked rows. No recommendation view.

CREATE OR REPLACE VIEW v_dual_evidence_genes AS
  SELECT * FROM read_csv_auto('results/tables/DUAL_EVIDENCE_CANDIDATES.csv') WHERE high_priority_research=TRUE;

CREATE OR REPLACE VIEW v_human_resting_zone_targets AS
  SELECT gene, height_direction, top_zone, gp_expression, zone_specificity, n_evidence_modes, source
  FROM read_csv_auto('results/tables/DUAL_EVIDENCE_CANDIDATES.csv')
  WHERE top_zone='Resting' AND gp_expression>0.4;

CREATE OR REPLACE VIEW v_fusion_capable_mature_length AS
  SELECT citation, species, fusion_capable_flag, endpoint, direction, source_location
  FROM studies WHERE is_final_length=TRUE AND fusion_capable_flag IN ('yes','mixed');

CREATE OR REPLACE VIEW v_direct_combinations AS
  SELECT intervention_a, intervention_b, relation FROM read_csv_auto('results/tables/PAIRWISE_INTERACTION_MATRIX.csv')
  WHERE relation <> 'untested/unknown';

CREATE OR REPLACE VIEW v_literature_adult_height AS
  SELECT pmid, title, year, journal, endpoint_type FROM literature
  WHERE endpoint_type IN ('adult_height','near_final') ORDER BY year DESC;

CREATE OR REPLACE VIEW v_zone_topology AS   -- any gene's zonal topology (genome-wide atlas)
  SELECT * FROM read_parquet('results/tables/GENOME_WIDE_STATE_SPECIFICITY.parquet');

-- unified text-search surface (studies + literature + interventions + calibration)
CREATE OR REPLACE VIEW v_search AS
  SELECT 'study' AS kind, citation AS title, notes AS text, source_location AS source FROM studies
  UNION ALL SELECT 'literature', title, endpoint_type, pmid FROM literature
  UNION ALL SELECT 'intervention', name, target||' | '||action, source_provenance FROM interventions
  UNION ALL SELECT 'calibration', intervention, adult_height_effect||' | translated='||surrogate_translated, source_id FROM surrogate_calibration
  UNION ALL SELECT 'access_failure', source, reason, request_route FROM access_failures;
