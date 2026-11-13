# WP P1C: CMUR calibration robustness

Completed 10 November 2026. This directory contains diagnostics only. No existing manuscript, data, CSV, figure, P1A result or Zenodo evidence was edited. No external data or literature was used. The uncertainty interpretation follows the Controller's fixed scientific basis.

## Authority and gates

Production authority is scripts/wp1b_r.py, whose OUT, load_panel, prepare, simulate, classify_array, short_metric and conventional_class functions locate and define the frozen inputs and mechanics. Its data/processed inputs are cmur_source_provenance_panel.csv, verified_anchor_map.csv, b2r_residual_pool.csv, b2r_short_window_support.csv and temporal_stability_comparison.csv. Medium-term rows come from the prior harmonized CSV, checked against the frozen branchA_medium_term.csv and branchB1_anchor_medium_term.csv located by production. The manuscript pipeline and scripts/verify_wp1a.py, scripts/verify_wp1b_r.py and manuscript/scripts/verify_integration.py establish production reproduction. No medium-term model was refitted.

Every mandatory P1A input was read: baseline_reproduction.csv, medium_term_harmonized_rules.csv, short_window_unresolved_causes.csv, short_window_variant_influence.csv, SCIENTIFIC_DIAGNOSTIC_SUMMARY.md, aggregate_summary.json, verification.json and independent_verification.json. Production gates reproduce 15/0/12 medium-term, 10/0/314 annual change, 2/0/295 turning, and 1/0/296 acceleration for Member States. P1A symmetric rules reproduce 12/0/15 and 15/0/12; both strict and exploratory midpoint removals reproduce zero new acceleration resolutions, for Member States and EU separately.

## Scenarios and paired simulation

- UNSCALED_STRESS (U) is PRIMARY: the empirical four-year pseudo-gap residual stress-test amplitude, with no claim to calibrated two-year errors or formal upper confidence bounds.
- BROWNIAN_SD_1_SQRT2 (B) is SECONDARY: exactly float(1/np.sqrt(2)) = 0.7071067811865475 multiplies each sampled residual before adding its midpoint and clipping. This is the Controller-specified Brownian-motion/bridge SD scenario, not the true CMUR correction.

The fixed seed is 20260926, with 20000 trajectories and the unchanged PCG64 SeedSequence([seed, geography index, variant index]) streams. Residual pool order, sampled identities, randomized signs, draw order, independent odd-year draws, clipping to [0,100], anchor values and country applicability are unchanged. A single identity/sign sequence feeds both scenarios; this also preserves draws with zero-valued or repeated residuals. Integer-index sampling was verified against production's value-sampling API, and every U path was array-equal to production simulate. The full regenerated U CSV is byte-identical to frozen production, not only its classes. Scenario B is evaluated only after this gate passes. No residual ×1/2 scenario was run.

paired_draw_audit.csv stores one row per geography/variant/eligible odd-year block, including original stream indices, hashes of identities, signs and unscaled/scaled pre-clipping residuals, and clipped-draw counts. Identical identity/sign hashes are deliberate: both scenarios use the same sampled objects. Verification independently replays the seeds, reconstructs these hashes and both paths, and checks B estimates and sign frequencies. No random stream is reseeded or rerouted after variant exclusion. Full trajectories are reproducible from the saved code and frozen inputs; they are not duplicated as a large draw archive.

The deterministic midpoint has one path, unchanged in B, and uses the production exact-zero acceleration representation for anchor–midpoint–anchor triples. It is structural information without invented probability values.

## Continuous support definitions

calibration_scenario_support.csv has 7,616 rows: two scenarios × 28 geographies × (12 annual + 11 turning + 11 acceleration) × four variants, through 2022 only. Applicable missing-coverage rows remain in the grid as unresolved; country_signed is NOT_APPLICABLE for CZ and EL. No endpoint claim is silently dropped from the production byte-identity check, but 2023–2024 are outside this WP's scenario/threshold analysis.

For stochastic annual change and acceleration, p_positive and p_negative are the empirical fractions strictly above/below zero; zeros support neither. preferred_sign/class selects the larger fraction, with TIE explicit. preferred_sign_support is max(p_positive,p_negative). distance_from_0_95_threshold is the signed difference support−0.95. Fractions describe the sensitivity ensemble; they are neither posterior probabilities nor calibrated frequentist confidence probabilities.

Turning is not a sign-valued metric: p_positive/p_negative and preferred_sign are intentionally blank. preferred_class is TURN/NO_TURN from the two constituent preferred signs (or TIE). Its preferred_sign_support is the MINIMUM of the two constituent max-sign supports, as indicated by support_definition. The earlier/later p_positive, p_negative, preferred sign, support and 0.95 class are exposed. At each threshold both constituent change classes are recomputed; unresolved either side means turning unresolved. source_aware_estimate retains the production strict reversal-event frequency separately. It is never relabelled as a turning sign-support probability. Deterministic rows have all probability/support quantities blank.

claim_level_continuous_support.csv has 1,904 scenario/claim rows. All applicable stochastic variants are considered, including applicable unresolved ones; NOT_APPLICABLE is omitted. Preferred-class status is MISSING_COVERAGE, PREFERRED_TIE, CONCORDANT or CONFLICT. Only CONCORDANT rows get weakest_variant_support = minimum preferred support across applicable variants. Conflicts receive no averaged probability. Turning concordance concerns TURN/NO_TURN, not equality of both underlying sign directions. This statistic does not assume variants are independent. The deterministic estimate/class is included solely as structural reference.

## Threshold and comparison conventions

The exact grid is 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 0.99. threshold_claim_classifications.csv contains 13,328 scenario/claim/threshold rows, with all four recomputed variant classes. Stochastic-only resolution requires at least TWO applicable stochastic variants, every one resolved and all the same class; normally three vote, but CZ and EL have two. A conflicting-direction claim has at least two incompatible resolved stochastic classes, even if another applicable vote is unresolved; conflict takes precedence in the exclusive status. All other claims are stochastic UNRESOLVED. Preferred conflicts below threshold remain distinct in the continuous output.

Production-style classification keeps the literal minimum of THREE applicable votes including deterministic_midpoint, with unanimous resolved classes across ALL applicable votes; any applicable unresolved vote vetoes resolution. Conventional classes remain frozen. Agreement/disagreement/unresolved is the production comparison with that conventional class. A source-aware resolution and a final comparison resolution are separate concepts. threshold_survival.csv has 84 rows and separately reports both classification families for Member States and EU27_2020. Conflicting stochastic status is distinct from a final disagreement with conventional direction.

scenario_classification_changes.csv retains all 952 geography/year/metric U–B pairs at 0.95, not only changes. It supplies both source-aware and final classes, new/lost resolution flags, preferred statuses and weakest supports where defined. Filter final_class_changed=True for every changed final identity. A paired support difference is supplied only where both scenarios have concordant preferred classes; it is descriptive and may compare different unanimous preferred labels if a flip were present.

continuous_support_distribution_summary.csv reports P10, quartiles, median, P90 and proportions ≥0.80/0.90/0.95/0.99. scenario_concordant uses each scenario's concordant subset and explicitly reports excluded conflict/tie/missing counts. common_concordant_pairs restricts support distributions to the same keys concordant in both scenarios, reducing selection-composition ambiguity. Quantiles use NumPy's default linear interpolation. Proportion denominators are the concordant subset, not the whole claim grid. Exclusion counts always describe the full grid. No significance test was performed.

## Exploratory multi-span design

multispan_pseudogap_diagnostic.csv covers all admissible t−2,t,t+2 four-year and t−4,t,t+4 eight-year triples within 2010–2022. Every central target is observed and admissible; no validation value was imputed. Four-year residual keys/values reproduce the frozen pool. The same-country/same-central-target paired subset requires the union of both triples (five admissible anchors). The Member-State exercise has broad country coverage adequate for descriptive summaries; the EU aggregate's three eight-year targets are displayed as sparse descriptive evidence only. Repeated/overlapping windows and countries are not independent observations.

multispan_summary.csv separates Member States and EU27_2020, all-eligible and paired subsets, with counts, geography coverage/list, median and mean |residual|, RMSE, raw signed-residual MAD and absolute-residual MAD. No normal-consistency multiplier is applied. multispan_paired_comparison.csv records exact-decimal absolute residual ordering for each pair from the released observation strings; this avoids floating-point false non-ties without adding a tolerance. Larger typical eight-year errors do not establish individual monotonicity or identify two-year errors. No exponent is fitted and no empirical two-year scale is extrapolated.

## Medium-term candidate and figures

medium_term_revision_candidate.csv contains four rule/population rows: symmetric two-sided 95% as candidate_PRIMARY and common directional 95% as SECONDARY_robustness. Counts, unresolved counts, median CI widths, paired width ratios and changed geographies are computed from verified P1A rows only. This does not approve a manuscript change.

threshold_survival.png/.pdf plots the whole-grid share lacking stochastic-only concordant resolution (solid) and final unresolved share (dotted); these overlap in the computed grid, demonstrating that the midpoint vote adds no veto there. continuous_support_distribution.png/.pdf shows ECDFs conditional on concordant preferred classes, with each denominator in the legend and missing/conflict exclusions documented here and in the CSV. Both are diagnostic candidates, not approved manuscript figures.

## Verification and reproduction

From the repository root:

```powershell
python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py baseline
python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py run
python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py report
```

The baseline stage runs the existing 44-check WP1A suite, 30-check WP1B-R suite (including full simulation rerun), and 17-check manuscript integration suite in an isolated disposable repository copy under this directory. All pass; every regenerated production CSV is byte-identical. Existing integration checks read the preserved LaTeX logs/PDF; no fresh manuscript build occurs. The run stage verifies U before B, produces/reloads all row data and independently verifies thresholds, aggregates, paired random replay, B sign frequencies and pseudo-gap admissibility. report regenerates every prose/table number from serialized outputs and creates figures, then repeats verification. aggregate can regenerate derived tables from the saved support CSV; verify can run checks alone.

input_integrity.json records all 371 pre-existing file hashes before and after, including the entire previous diagnostic directory; zero changes and zero additions outside this new directory are required. Verification never writes to a pre-existing path. The supplied workspace has no .git metadata, so Git diff/commit is unavailable; SHA-256 provides whole-file integrity evidence. All deliverables are saved locally. No blocker remains.

## Exact file inventory

- README.md
- SCIENTIFIC_DECISION_SUMMARY.md
- baseline_gate.csv
- calibration_scenario_support.csv
- claim_level_continuous_support.csv
- continuous_support_distribution.pdf
- continuous_support_distribution.png
- continuous_support_distribution_summary.csv
- input_integrity.json
- medium_term_revision_candidate.csv
- multispan_paired_comparison.csv
- multispan_pseudogap_diagnostic.csv
- multispan_summary.csv
- paired_draw_audit.csv
- run_diagnostics.py
- scenario_classification_changes.csv
- scenario_u_production_reproduction.csv
- threshold_claim_classifications.csv
- threshold_survival.csv
- threshold_survival.pdf
- threshold_survival.png
- verification.json
