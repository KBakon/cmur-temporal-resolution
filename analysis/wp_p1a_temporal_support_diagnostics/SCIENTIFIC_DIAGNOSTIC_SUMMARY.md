# Scientific diagnostic summary

Diagnostics of the frozen CMUR analysis, 10 November 2026. No manuscript claim or authoritative result was changed. A / D / U means agreement / resolved-class disagreement / unresolved, comparing the conventional and source-aware classes. Headline counts cover the fixed Member-State grid through 2022; EU27_2020 and endpoints are reported separately.

1. **The medium-term absence of conflicting resolved direction survives both identical inference rules.** EU27_2020 remains positive in A and B1 under both rules.

| Rule | Member States: agreement / disagreement / unresolved | EU27_2020 |
| --- | --- | --- |
| COMMON_TWO_SIDED_95 | 12 / 0 / 15 | 1 / 0 / 0 |
| COMMON_DIRECTIONAL_95 | 15 / 0 / 12 | 1 / 0 / 0 |

Under COMMON_TWO_SIDED_95, B1 changes from positive to unresolved for Belgium (BE), and from negative to unresolved for Finland (FI) and Luxembourg (LU); all three comparisons change from agreement to unresolved. A does not change. Under COMMON_DIRECTIONAL_95, A changes from unresolved to positive for Slovenia (SI); B1 remains unresolved there, so no comparison changes. These are all changed geographies.

2. **The asymmetric rules partly conceal the loss of directional precision.** The current A/B1 unresolved counts are 8/12, a gap of four. Common two-sided inference gives 8/15, a gap of seven: three of those seven additional B1 unresolved classifications (42.9%) are masked by B1's more permissive current directional rule. Common directional inference gives 7/12, a gap of five. Thus the asymmetry also understates that same-rule gap by one. These are classification differences, not a causal decomposition of sampling uncertainty. B1's standard 95% interval is wider for all 27 Member States; median widths are 0.3386 (A) and 0.6219 (B1) percentage points/year, with a median paired width ratio of 1.858. Harmonization changes classes, not slopes or intervals.

3. **Most current short-window unresolved claims have estimable stochastic variants that fail their sign-support rules.** Coverage-first exclusive attribution reconciles exactly:

| Metric | Unresolved | Coverage primary | Stochastic primary | Estimable midpoint unresolved (overlap) | Structural zero (overlap) |
| --- | --- | --- | --- | --- | --- |
| annual_change | 314 | 8 (2.5%) | 306 (97.5%) | 8 (2.5%) | 0 |
| turning | 295 | 9 (3.1%) | 286 (96.9%) | 8 (2.7%) | 0 |
| acceleration | 296 | 9 (3.0%) | 287 (97.0%) | 161 (54.4%) | 158 |

Every through-2022 unresolved Member-State claim has at least one applicable stochastic variant unresolved, including coverage failures. Once coverage failures are separated, the stochastic primary counts are 306, 286 and 287. The midpoint is additionally unresolved in 16, 17 and 170 claims respectively when missing-coverage rows are included; the estimable counts in the table isolate geometry from absent estimates. Structural zero accounts for 158 acceleration claims, or 53.4% of unresolved acceleration claims. All midpoint failures overlap stochastic failures; none is a sole veto that prevents otherwise unanimous stochastic resolution. Fewer-than-three applicability and consensus-only failure each occur zero times. One acceleration claim (CY, 2017) contains conflicting resolved variant signs (country negative, midpoint positive), alongside unresolved pooled stochastic variants; it is not a resolved disagreement between final branches.

EU27_2020 has 12/11/11 unresolved annual-change/turning/acceleration claims: all have estimable stochastic uncertainty and none lacks coverage. Six acceleration claims also have midpoint structural zeros (54.5%); other midpoint-unresolved counts are zero. For 2023–2024, every metric has 54 Member-State and two EU endpoint claims, all unresolved because of absent required coverage. These endpoints are excluded from headline denominators. The flags overlap and cannot be added as independent causal proportions. The exclusive hierarchy allocates overlaps to the first applicable cause; its zero midpoint-primary count does not mean midpoint geometry is absent.

4. **The qualitative short-window conclusion survives every leave-one-variant-out diagnostic.** Counts are:

| Voting diagnostic | Annual change A / D / U | Turning A / D / U | Acceleration A / D / U |
| --- | --- | --- | --- |
| LOO_pooled_signed | 10 / 0 / 314 | 2 / 0 / 295 | 1 / 0 / 296 |
| LOO_country_signed | 20 / 0 / 304 | 12 / 0 / 285 | 5 / 0 / 292 |
| LOO_pooled_absolute_random_sign | 10 / 0 / 314 | 2 / 0 / 295 | 1 / 0 / 296 |
| LOO_deterministic_midpoint | 10 / 0 / 314 | 2 / 0 / 295 | 1 / 0 / 296 |
| EXPLORATORY_midpoint_structural_reference | 10 / 0 / 314 | 2 / 0 / 295 | 1 / 0 / 296 |

Only removing country_signed changes claims: 10 additional annual-change agreements (3.2% of the current 314 unresolved), 10 turning agreements (3.4% of 295), and four acceleration agreements (1.4% of 296). All 24 changed comparisons move from unresolved to agreement; their exact geography/year/metric keys and classes are in short_window_variant_influence.csv. Even this diagnostic leaves 93.8%, 96.0% and 98.3% of the respective Member-State grids unresolved. Removing either pooled variant or the midpoint changes no source-aware or final comparison class. Removing the midpoint resolves **zero additional acceleration claims**, under both the strict minimum-three LOO and the explicitly EXPLORATORY stochastic-only minimum-two convention. EU27_2020 remains 0/0/12, 0/0/11 and 0/0/11 under every diagnostic. Endpoints remain 0/0/54 for Member States and 0/0/2 for EU27_2020 for each metric and diagnostic.

5. **The temporal-resolution support boundary is genuinely mixed within the tested design.** Missing or flagged anchors and absent endpoint coverage create direct provenance constraints. Within covered windows, the four-year calibration residuals, fixed stochastic sign-support thresholds, and unanimity requirement jointly determine support. Country-variant removal shows a modest consensus sensitivity, while the midpoint's frequent geometric zeros do not independently account for the observed acceleration unresolvedness. The unchanged majority-unresolved outcome across these diagnostics supports the qualitative boundary under this calibrated design. These diagnostics hold calibration and sign-support thresholds fixed; they do not identify a separate causal percentage for calibration versus stochastic threshold versus consensus, or establish the uncertainty distribution of true annual CMUR errors.

The single most decision-relevant finding is that removing the deterministic midpoint resolves no additional acceleration claims: stochastic support remains limiting, while the medium-term zero-conflict finding survives harmonized inference.
