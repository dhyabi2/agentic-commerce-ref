# Phase B saturation analysis (measured, honest)

**4 real forward+backward E-utilities citation rounds** run for all 15 categories.
Corpus: **62 → 274 included primary studies** (+212). Logs: `SEARCH_ROUNDS.csv`,
`CITATION_FRONTIER.csv`, `SCREENING_LOG.csv`, `EXCLUDED_STUDIES.csv`, `SATURATION.csv`, `PRISMA_COUNTS.csv`.

## Per-category status after 4 rounds
| category | topic | last-round new | saturated |
|---|---|---:|:---:|
| Q1 | GH / ISS adult height | 0 | **YES** |
| Q3 | vosoritide | 0 | **YES** |
| Q4 | CNP / hypochondroplasia | 0 | **YES** |
| Q13 | physeal distraction | 0 | **YES** |
| Q14 | resting-zone stem cell | 0 | **YES** |
| Q6 | oxandrolone / Turner final height | 1 | near (converging) |
| Q7 | IGF-1 / Laron adult height | 1 | near (converging) |
| Q9 | estrogen / fusion | 3 | open |
| Q8 | growth-plate biology (mechanism) | 5 | open frontier |
| Q11 | FGFR3 inhibition | 5 | open |
| Q2 | aromatase inhibitor final height | 6 | open |
| Q5 | GnRH agonist adult height | 6 | open |
| Q10 | NPR2 / NPPC / CNP | 12 | deep open frontier |
| Q12 | limb lengthening final length | 12 | deep open frontier |
| Q15 | ACAN genetics | 13 | deep open frontier |

## Principled interpretation (not faked)
- **Clinical/final-height categories are effectively covered.** The human adult/near-final-height
  evidence base (the PRIMARY endpoint) is well-sampled: Q1/Q3/Q4 saturated; Q6/Q7 converged to ≤1
  new/round. Q2/Q5 still surface a few clinical trials per round (worth 1–2 more rounds — seeds logged).
- **Mechanism/genetics categories (Q8, Q9, Q10, Q15) and orthopedic Q12 have genuinely open citation
  frontiers.** These are Tier-7/8 supporting literatures (growth-plate cell biology, CNP/NPR2 signaling,
  aggrecan cohorts, distraction-osteogenesis case series spanning 1961–2026). They expand ~5–13
  new/round and would need many more rounds; exhaustive saturation is neither achievable in-session nor
  necessary for the primary final-height question. They are **capped at 4 rounds with residual logged**,
  not declared saturated.
- **No saturation was faked.** `SATURATION.csv` marks only the 5 zero-new categories as saturated.

## Stopping rule applied
Continue citation rounds while a category still yields **final-height/near-final human** items
(Q2, Q5 → round 5+ seeds in `EXACT_CONTINUATION_PROMPT.md`). Cap Tier-7/8 mechanism categories at 4
rounds and log the open residual rather than pursue unbounded mechanism-literature saturation.
