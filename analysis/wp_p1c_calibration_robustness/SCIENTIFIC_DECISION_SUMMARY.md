# Scientific decision evidence: calibration robustness

10 November 2026. U is the PRIMARY empirical unscaled four-year stress test. B is the SECONDARY Brownian-model amplitude scenario, multiplying the same sampled residuals by 1/√2. Neither scenario supplies calibrated two-year uncertainty probabilities. A/D/U means agreement / resolved-class disagreement / unresolved against the released series. This report records diagnostics; it makes no manuscript, publication, novelty, P4 or journal decision.

## 1. Does the short-window boundary survive Brownian scaling?

Yes, within these two scenarios and the pre-specified threshold grid. Most claims lack resolved cross-variant sign support even after amplitude reduction. The statement concerns CMUR under the tested calibration families and consensus design, not the unknown distribution of true annual errors.

## 2. What changes at the 0.95 reference threshold?

| Population | Metric | U A/D/U | B A/D/U | New resolutions | Lost | Change unresolved, pp |
| --- | --- | --- | --- | --- | --- | --- |
| Member States | annual_change | 10/0/314 | 29/0/295 | 19 | 0 | -5.86 |
| Member States | turning | 2/0/295 | 13/0/284 | 11 | 0 | -3.70 |
| Member States | acceleration | 1/0/296 | 7/0/290 | 6 | 0 | -2.02 |
| EU27_2020 | annual_change | 0/0/12 | 0/0/12 | 0 | 0 | +0.00 |
| EU27_2020 | turning | 0/0/11 | 0/0/11 | 0 | 0 | +0.00 |
| EU27_2020 | acceleration | 0/0/11 | 0/0/11 | 0 | 0 | +0.00 |

All 36 changed final comparisons become agreements; 0 lose resolution. Exact geography/year/metric identities and old/new classes are in scenario_classification_changes.csv (filter final_class_changed=True). Scaling changes a meaningful set of individual classifications but reduces unresolvedness by only the percentage-point amounts above. Annual change is most sensitive in absolute counts and percentage-point change; acceleration is least sensitive by those measures. Relative growth from tiny baseline agreement counts gives a different ranking and should not be confused with panel-level sensitivity. EU27_2020 gains no resolution.

## 3. Is the conclusion dominated by choosing 0.95?

No majority conclusion depends uniquely on 0.95: every metric remains more than half unresolved throughout 0.80, 0.85, 0.90, 0.925, 0.95, 0.975 and 0.99, in both scenarios. Detailed counts change with threshold:

| Scenario | Metric | Final unresolved range over grid | Stochastic concordant at .80 | At .99 |
| --- | --- | --- | --- | --- |
| U: unscaled stress | annual_change | 85.2–99.1% | 48 | 3 |
| U: unscaled stress | turning | 93.6–99.7% | 19 | 1 |
| U: unscaled stress | acceleration | 93.6–100.0% | 19 | 0 |
| B: Brownian SD | annual_change | 75.6–96.3% | 79 | 12 |
| B: Brownian SD | turning | 86.5–98.0% | 40 | 6 |
| B: Brownian SD | acceleration | 89.6–100.0% | 31 | 0 |

All threshold-level stochastic conflict counts are 0: the opposing preferred classes described below do not both reach these thresholds. Stochastic-only concordance and production-style final resolution counts coincide in these simulations across this grid; keeping the midpoint vote adds no veto at these thresholds. The grid is a robustness diagnosis, not a threshold-selection exercise. Resolved disagreement does occur under B at low thresholds: BE 2019 annual_change, thresholds 0.8, 0.85. This annual change is conventional negative versus source-aware positive; there is no final resolved disagreement at 0.95. EU27_2020 remains entirely unresolved throughout the grid in both scenarios.

## 4–5. What do continuous supports reveal?

The weakest-variant statistic is the minimum preferred-sign support across applicable stochastic variants only when their preferred classes agree. It is a robustness statistic, not a joint or calibrated probability. For turning it is the minimum of constituent-change supports across variants, with the preferred TURN/NO_TURN class inferred from the constituent preferred signs. Strict reversal-event frequencies are stored separately and never substituted for this support.

Member-State distributions conditional on concordant preferred classes are:

| Scenario | Metric | n concordant | Preferred conflicts | P10 | Q25 | Median | Q75 | P90 | ≥.80 | ≥.90 | ≥.95 | ≥.99 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| U: unscaled stress | annual_change | 263 | 53 | 0.537 | 0.597 | 0.632 | 0.797 | 0.805 | 18.3% | 6.1% | 3.8% | 1.1% |
| U: unscaled stress | turning | 209 | 79 | 0.524 | 0.590 | 0.602 | 0.736 | 0.799 | 9.1% | 1.9% | 1.0% | 0.5% |
| U: unscaled stress | acceleration | 152 | 136 | 0.479 | 0.515 | 0.557 | 0.666 | 0.801 | 12.5% | 0.7% | 0.7% | 0.0% |
| B: Brownian SD | annual_change | 276 | 40 | 0.563 | 0.601 | 0.709 | 0.801 | 0.956 | 28.6% | 11.2% | 10.5% | 4.3% |
| B: Brownian SD | turning | 226 | 62 | 0.549 | 0.598 | 0.649 | 0.797 | 0.834 | 17.7% | 5.8% | 5.8% | 2.7% |
| B: Brownian SD | acceleration | 161 | 127 | 0.479 | 0.518 | 0.567 | 0.697 | 0.882 | 19.3% | 7.5% | 4.3% | 0.0% |

The distribution denominators differ because scaling can change preferred-class concordance. Missing-coverage counts are unchanged (8/9/9 by metric); no preferred ties occur. Conflicting preferred classes are explicitly excluded from a single weakest-support number, not averaged. The common-concordant paired subset is also saved: its medians U → B are annual_change: 0.632 → 0.722; turning: 0.602 → 0.669; acceleration: 0.557 → 0.574. EU27_2020 conditional medians U → B are annual_change: 0.550 → 0.612; turning: 0.527 → 0.547; acceleration: 0.527 → 0.553; none reaches cross-variant support of 0.80. These graded shifts and numerous weak direction conflicts are obscured by a single 0.95 class. Acceleration shows especially low support and frequent preferred-direction conflicts.

## 6. What does the exploratory eight-year exercise show?

Task F is feasible for a descriptive Member-State exercise. All predictions use genuinely admissible outer and target anchors. Values below are in CMUR percentage points; MAD is the unscaled median absolute deviation about the median signed residual.

| Population | Subset | Span | n | Geographies | Median \|e\| | Mean \|e\| | RMSE | Signed MAD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Member States | all_eligible | 4 | 130 | 27 | 0.600 | 1.004 | 1.563 | 0.550 |
| Member States | all_eligible | 8 | 77 | 27 | 0.750 | 1.432 | 2.085 | 0.850 |
| Member States | paired_same_geo_centre | 4 | 76 | 26 | 0.525 | 0.918 | 1.333 | 0.500 |
| Member States | paired_same_geo_centre | 8 | 76 | 26 | 0.750 | 1.432 | 2.092 | 0.825 |
| EU27_2020 | all_eligible | 4 | 5 | 1 | 0.100 | 0.170 | 0.201 | 0.150 |
| EU27_2020 | all_eligible | 8 | 3 | 1 | 0.300 | 0.233 | 0.268 | 0.050 |
| EU27_2020 | paired_same_geo_centre | 4 | 3 | 1 | 0.100 | 0.150 | 0.185 | 0.150 |
| EU27_2020 | paired_same_geo_centre | 8 | 3 | 1 | 0.300 | 0.233 | 0.268 | 0.050 |

The paired Member-State subset has the same geography and central target for both spans, with all five required anchors admissible. Eight-year absolute residuals are larger in 54 cases, smaller in 20, and equal in 2 (ordering evaluated at exact released decimal precision). Thus pooled and paired aggregate scales are larger for eight years, consistent with smaller typical error at shorter observed spans, but casewise monotonicity fails. Nested windows, overlapping targets and within-country dependence prevent treating these as independent replications. EU27_2020 has only 3 pairs and is descriptive only. This comparison does not observe two-year validation errors, validate Brownian dynamics or the 1/√2 factor, calibrate a two-year scale, or support a validated scaling exponent.

## 7. What level of scientific interpretation is supported?

The evidence goes beyond a threshold-specific result: the majority boundary survives the full threshold grid and the approved model-motivated amplitude reduction, with continuous summaries showing predominantly weak support. It supports a **robust CMUR-specific temporal-support boundary within this tested sensitivity design**. The computational exercise also illustrates a method for separating coverage, preferred-direction conflicts, support thresholds and amplitude assumptions; it does not establish a general empirical boundary for other indicators. Calibration from four-year residuals to actual two-year error remains unidentified without further assumptions or validation. The Controller retains all final contribution, novelty and P4 decisions.

## Medium-term revision candidate (no new models)

| Candidate role | Population | A/D/U | A/B1 unresolved | Median A/B1 CI width | Median paired ratio | Changed geographies |
| --- | --- | --- | --- | --- | --- | --- |
| candidate_PRIMARY | Member States | 12/0/15 | 8/15 | 0.3386/0.6219 | 1.858 | A: none; B1: BE;FI;LU |
| candidate_PRIMARY | EU27_2020 | 1/0/0 | 0/0 | 0.0975/0.1016 | 1.042 | A: none; B1: none |
| SECONDARY_robustness | Member States | 15/0/12 | 7/12 | 0.3386/0.6219 | 1.858 | A: SI; B1: none |
| SECONDARY_robustness | EU27_2020 | 1/0/0 | 0/0 | 0.0975/0.1016 | 1.042 | A: none; B1: none |

COMMON_TWO_SIDED_95 is the candidate symmetric primary; COMMON_DIRECTIONAL_95 is secondary robustness. BE changes B1 positive → unresolved; FI and LU negative → unresolved under the candidate primary, changing their comparisons to unresolved. SI changes A unresolved → positive under the secondary, while its comparison remains unresolved. These rows are reused from verified P1A fits; no medium-term model was refitted here.

## Verification

Production and P1A gates passed before alternatives. Scenario U's complete production support CSV is byte-identical (including endpoints and all stochastic fields/classes); B alone changes residual amplitude. Independent seed replay verifies residual identities/signs/order and every stored B sign support. All 13328 threshold claim rows are independently reconstructed from serialized continuous support, and all reported aggregate tables are generated from serialized files. Existing verification suites passed in isolation; 371 pre-existing files, including all P1A outputs, manuscript, figures, frozen sources and Zenodo evidence, remain byte-identical. No scientific or technical blocker remains; the workspace lacks Git metadata, so integrity is verified by SHA-256 rather than a Git diff. Existing integration checks use the preserved build logs/PDF, not a new LaTeX build.

The most decision-relevant finding is that Brownian scaling raises 0.95 resolutions but leaves the majority short-window boundary intact even at 0.80, so the qualitative finding is not an artifact of the 0.95 cutoff alone.
