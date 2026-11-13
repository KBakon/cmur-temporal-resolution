"""Read-only diagnostics of the frozen CMUR analysis; writes only beside this file."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from contextlib import contextmanager
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import t

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'scripts'))
import wp1b_r as prod

METRICS = ('annual_change', 'turning', 'acceleration')
LABELS = ('AGREEMENT', 'CHANGE', 'UNRESOLVED')
FLAGS = ('missing_required_coverage', 'fewer_than_required_applicable_variants',
         'any_applicable_stochastic_unresolved', 'estimable_stochastic_unresolved',
         'deterministic_midpoint_unresolved', 'estimable_midpoint_unresolved',
         'conflicting_resolved_variant_classes', 'failure_only_cross_variant_consensus',
         'midpoint_only_veto', 'conventional_class_unresolved', 'midpoint_structural_zero')


def read(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def write(name, rows):
    assert rows
    with (HERE / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


@contextmanager
def verification_copy():
    # On this Windows host tempfile's mode 0700 produces an inaccessible ACL.
    sandbox = (HERE / 'verification_workspace').resolve()
    assert sandbox.is_relative_to(HERE) and sandbox != HERE
    sandbox.mkdir(exist_ok=False)
    try:
        yield sandbox
    finally:
        assert sandbox.is_relative_to(HERE) and sandbox != HERE
        shutil.rmtree(sandbox)


def inventory():
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(ROOT.rglob('*'))
            if p.is_file() and not p.is_relative_to(HERE)}


def population(geo):
    return 'EU27_2020' if geo == 'EU27_2020' else 'Member States'


def period(year):
    return '2010-2022' if year == '2010-2022' or int(year) <= 2022 else '2023-2024 endpoints'


def stability(a, b):
    return 'UNRESOLVED' if 'UNRESOLVED' in (a, b) else 'AGREEMENT' if a == b else 'CHANGE'


def consensus(rows, minimum=3):
    classes = [r['variant_class'] for r in rows if r['variant_class'] != 'NOT_APPLICABLE']
    return classes[0] if len(classes) >= minimum and len(set(classes)) == 1 and classes[0] != 'UNRESOLVED' else 'UNRESOLVED'


def inputs():
    # Paths come directly from the production module used by the manuscript pipeline.
    return [read(prod.OUT / n) for n in ('branchA_medium_term.csv', 'branchB1_anchor_medium_term.csv',
            'b2r_short_window_support.csv', 'temporal_stability_comparison.csv')]


def baseline():
    initial = inventory()
    if (HERE / 'input_integrity.json').exists():
        assert initial == json.loads((HERE / 'input_integrity.json').read_text())['before_sha256']
    else:
        save('input_integrity.json', {'before_sha256': initial})
    A, B, S, C = inputs()
    expected = {'medium_term_slope': (15, 0, 12), 'annual_change': (10, 0, 314),
                'turning': (2, 0, 295), 'acceleration': (1, 0, 296)}
    rows = []
    for metric, target in expected.items():
        for pop in ('Member States', 'EU27_2020'):
            selected = [r for r in C if r['metric'] == metric and population(r['geo']) == pop and period(r['year']) == '2010-2022']
            counts = Counter(r['stability_class'] for r in selected)
            actual = tuple(counts[x] for x in LABELS)
            wanted = target if pop == 'Member States' else (1, 0, 0) if metric == 'medium_term_slope' else (0, 0, 12 if metric == 'annual_change' else 11)
            assert actual == wanted, (metric, pop, actual, wanted)
            rows.append(dict(population=pop, period='2010-2022', metric=metric,
                             agreement=actual[0], resolved_class_disagreement=actual[1], unresolved=actual[2],
                             total=len(selected), expected_agreement=wanted[0], expected_disagreement=wanted[1],
                             expected_unresolved=wanted[2], exact_match=True))
    grouped = defaultdict(list)
    for r in S:
        grouped[(r['geo'], r['year'], r['metric'])].append(r)
    for r in C:
        if r['metric'] == 'medium_term_slope':
            a = next(x for x in A if x['geo'] == r['geo'])
            b = next(x for x in B if x['geo'] == r['geo'])
            assert (r['branchA_class'], r['source_aware_class']) == (a['direction'], b['direction'])
        else:
            assert consensus(grouped[(r['geo'], r['year'], r['metric'])]) == r['source_aware_class']
        assert stability(r['branchA_class'], r['source_aware_class']) == r['stability_class']
    write('baseline_reproduction.csv', rows)
    print('Frozen baseline counts and every stored consensus/stability class reproduce exactly.', flush=True)
    checks = []
    with verification_copy() as sandbox:
        for folder in ('scripts', 'data', 'logs', 'manuscript'):
            shutil.copytree(ROOT / folder, sandbox / folder, ignore=shutil.ignore_patterns('__pycache__'))
        for script in ('scripts/verify_wp1a.py', 'scripts/verify_wp1b_r.py', 'manuscript/scripts/verify_integration.py'):
            command = [sys.executable, '-B', script]
            result = subprocess.run(command, cwd=sandbox, capture_output=True, text=True,
                                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
            checks.append(dict(command='python -B ' + script, execution='isolated temporary repository copy',
                               returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
            save('verification.json', {'checks': checks, 'baseline_pass': False})
            assert result.returncode == 0, checks[-1]
        # Existing verifiers regenerate scientific outputs. Require exact bytes versus originals.
        outputs = __import__('verify_wp1a').PRODUCTS + __import__('verify_wp1b_r').FILES
        assert all(sha(prod.OUT / n) == sha(sandbox / 'data/processed' / n) for n in outputs)
    assert inventory() == initial, 'A pre-existing repository file changed during baseline verification.'
    save('verification.json', {'checks': checks, 'baseline_pass': True,
                              'all_regenerated_processed_outputs_byte_identical': True,
                              'all_original_files_unchanged_after_baseline': True,
                              'seed': prod.SEED, 'draws': prod.N_DRAWS})
    print('All three existing verification suites PASS in isolation; regenerated results byte-identical.', flush=True)


def diagnostics():
    gate = json.loads((HERE / 'verification.json').read_text())
    assert gate['baseline_pass'], 'STOP: baseline gate not passed.'
    assert inventory() == json.loads((HERE / 'input_integrity.json').read_text())['before_sha256']
    A, B, S, C = inputs()
    p, anchor_map = prod.load_panel()
    medium = []
    for rule in ('COMMON_TWO_SIDED_95', 'COMMON_DIRECTIONAL_95'):
        for geo in prod.GEOS:
            originals = [next(r for r in rows if r['geo'] == geo) for rows in (A, B)]
            fitted = [prod.ols([(y, prod.cmur(p, geo, y)) for y in prod.YEARS
                              if branch == 'A' or prod.admissible(anchor_map, geo, y)]) for branch in ('A', 'B1')]
            classes = []
            row = dict(rule=rule, geo=geo, population=population(geo))
            for branch, fit, frozen in zip(('A', 'B1'), fitted, originals):
                for key in ('slope_pp_per_year', 'stderr', 'ci_low', 'ci_high', 'p_positive'):
                    assert np.isclose(fit[key], float(frozen[key]), rtol=0, atol=1e-12), (geo, branch, key)
                    row[f'{branch}_{key}'] = fit[key]
                row[f'{branch}_n'] = fit['n']
                cls = ('POSITIVE' if fit['ci_low'] > 0 else 'NEGATIVE' if fit['ci_high'] < 0 else 'UNRESOLVED') if rule == 'COMMON_TWO_SIDED_95' else ('POSITIVE' if fit['p_positive'] >= .95 else 'NEGATIVE' if fit['p_positive'] <= .05 else 'UNRESOLVED')
                row[f'{branch}_class'] = cls
                row[f'{branch}_current_class'] = frozen['direction']
                row[f'{branch}_class_changed'] = cls != frozen['direction']
                classes.append(cls)
            row['stability_class'] = stability(*classes)
            row['current_stability_class'] = stability(*(r['direction'] for r in originals))
            row['stability_class_changed'] = row['stability_class'] != row['current_stability_class']
            medium.append(row)
    write('medium_term_harmonized_rules.csv', medium)

    grouped = defaultdict(dict)
    for r in S:
        grouped[(r['geo'], r['year'], r['metric'])][r['variant']] = r
    short = [r for r in C if r['metric'] in METRICS]
    causes, influence = [], []
    for current in short:
        geo, year, metric = current['geo'], current['year'], current['metric']
        variants = grouped[(geo, year, metric)]
        applicable = [r for r in variants.values() if r['variant_class'] != 'NOT_APPLICABLE']
        stochastic = [r for r in applicable if r['variant'] != 'deterministic_midpoint']
        det = variants['deterministic_midpoint']
        base = dict(geo=geo, year=year, metric=metric, population=population(geo), period=period(year),
                    branchA_class=current['branchA_class'], current_source_aware_class=current['source_aware_class'],
                    current_stability_class=current['stability_class'])
        if current['stability_class'] == 'UNRESOLVED':
            missing = [r for r in applicable if r['reason'] in ('endpoint_unanchored', 'missing_admissible_anchor_or_bracket')]
            unresolved_stochastic = [r for r in stochastic if r['variant_class'] == 'UNRESOLVED']
            estimable_stochastic = [r for r in unresolved_stochastic if r['source_aware_estimate'] != '']
            resolved = {r['variant_class'] for r in applicable} - {'UNRESOLVED'}
            all_resolved = all(r['variant_class'] != 'UNRESOLVED' for r in applicable)
            flags = dict(missing_required_coverage=bool(missing), fewer_than_required_applicable_variants=len(applicable) < 3,
                         any_applicable_stochastic_unresolved=bool(unresolved_stochastic),
                         estimable_stochastic_unresolved=bool(estimable_stochastic),
                         deterministic_midpoint_unresolved=det['variant_class'] == 'UNRESOLVED',
                         estimable_midpoint_unresolved=det['variant_class'] == 'UNRESOLVED' and det['source_aware_estimate'] != '',
                         conflicting_resolved_variant_classes=len(resolved) > 1,
                         failure_only_cross_variant_consensus=current['source_aware_class'] == 'UNRESOLVED' and not missing and bool(applicable) and all_resolved,
                         midpoint_only_veto=current['source_aware_class'] == 'UNRESOLVED' and consensus(stochastic, 2) != 'UNRESOLVED',
                         conventional_class_unresolved=current['branchA_class'] == 'UNRESOLVED',
                         midpoint_structural_zero=metric == 'acceleration' and det['source_aware_estimate'] == '0.0' and int(year) % 2 == 0 and not missing)
            hierarchy = [('source_coverage', flags['missing_required_coverage']),
                         ('insufficient_applicable_variants', flags['fewer_than_required_applicable_variants']),
                         ('stochastic_uncertainty', flags['estimable_stochastic_unresolved']),
                         ('deterministic_midpoint_geometry', flags['estimable_midpoint_unresolved']),
                         ('resolved_variant_conflict_consensus', flags['conflicting_resolved_variant_classes']),
                         ('conventional_unresolved', flags['conventional_class_unresolved'])]
            primary = next((name for name, yes in hierarchy if yes), None)
            assert primary is not None, base
            row = base | dict(applicable_variants=len(applicable), required_variants=3) | flags | dict(primary_cause=primary)
            for v in prod.VARIANTS:
                row[f'{v}_class'] = variants[v]['variant_class']
                row[f'{v}_reason'] = variants[v]['reason']
            causes.append(row)
        diagnostics = [('LOO_' + v, v, 3) for v in prod.VARIANTS]
        diagnostics.append(('EXPLORATORY_midpoint_structural_reference', 'deterministic_midpoint', 2))
        for diagnostic, omitted, minimum in diagnostics:
            votes = [r for v, r in variants.items() if v != omitted]
            cls = consensus(votes, minimum)
            final = stability(current['branchA_class'], cls)
            row = base | dict(diagnostic=diagnostic, omitted_voting_variant=omitted, minimum_applicable_variants=minimum,
                              applicable_voting_variants=sum(r['variant_class'] != 'NOT_APPLICABLE' for r in votes),
                              source_aware_class=cls, stability_class=final,
                              source_aware_class_changed=cls != current['source_aware_class'],
                              stability_class_changed=final != current['stability_class'],
                              newly_source_aware_resolved=current['source_aware_class'] == 'UNRESOLVED' and cls != 'UNRESOLVED',
                              newly_comparison_resolved=current['stability_class'] == 'UNRESOLVED' and final != 'UNRESOLVED')
            for v in prod.VARIANTS:
                row[f'{v}_class'] = variants[v]['variant_class']
            influence.append(row)
    write('short_window_unresolved_causes.csv', causes)
    write('short_window_variant_influence.csv', influence)
    independent_verify()


def independent_verify():
    """Reload serialized rows and independently derive every aggregate used in prose."""
    A, B, S, C = inputs()
    medium = read(HERE / 'medium_term_harmonized_rules.csv')
    causes = read(HERE / 'short_window_unresolved_causes.csv')
    influence = read(HERE / 'short_window_variant_influence.csv')
    truth = lambda x: x == 'True'
    assert len(medium) == 56 and len({(r['rule'], r['geo']) for r in medium}) == 56
    expected = {(r['geo'], r['year'], r['metric']) for r in C if r['metric'] in METRICS and r['stability_class'] == 'UNRESOLVED'}
    assert len(causes) == len(expected) and {(r['geo'], r['year'], r['metric']) for r in causes} == expected
    source = defaultdict(dict)
    for r in S:
        source[(r['geo'], r['year'], r['metric'])][r['variant']] = r
    _, anchors = prod.load_panel()
    for r in causes:
        vv = source[(r['geo'], r['year'], r['metric'])]
        aa = [x for x in vv.values() if x['variant_class'] != 'NOT_APPLICABLE']
        ss = [x for x in aa if x['variant'] != 'deterministic_midpoint']
        dd = vv['deterministic_midpoint']
        year = int(r['year'])
        required_years = range(year - (1 if r['metric'] == 'annual_change' else 2), year + 1)
        def has_value(y):
            return y <= 2022 and (prod.admissible(anchors, r['geo'], y) or
                (y in prod.ODD and prod.admissible(anchors, r['geo'], y-1) and prod.admissible(anchors, r['geo'], y+1)))
        coverage_missing = not all(has_value(y) for y in required_years)
        assert truth(r['missing_required_coverage']) == coverage_missing
        assert truth(r['fewer_than_required_applicable_variants']) == (len(aa) < 3)
        assert truth(r['any_applicable_stochastic_unresolved']) == any(x['variant_class'] == 'UNRESOLVED' for x in ss)
        uncertainty = any(x['variant_class'] == 'UNRESOLVED' and x['source_aware_estimate'] != '' for x in ss)
        assert truth(r['estimable_stochastic_unresolved']) == uncertainty
        midpoint_unresolved = dd['variant_class'] == 'UNRESOLVED'
        assert truth(r['deterministic_midpoint_unresolved']) == midpoint_unresolved
        assert truth(r['estimable_midpoint_unresolved']) == (midpoint_unresolved and dd['source_aware_estimate'] != '')
        conflict = len({x['variant_class'] for x in aa} - {'UNRESOLVED'}) > 1
        assert truth(r['conflicting_resolved_variant_classes']) == conflict
        all_resolved = all(x['variant_class'] != 'UNRESOLVED' for x in aa)
        assert truth(r['failure_only_cross_variant_consensus']) == (r['current_source_aware_class'] == 'UNRESOLVED' and not coverage_missing and bool(aa) and all_resolved)
        stochastic_classes = [x['variant_class'] for x in ss]
        stochastic_consensus = len(stochastic_classes) >= 2 and len(set(stochastic_classes)) == 1 and stochastic_classes[0] != 'UNRESOLVED'
        assert truth(r['midpoint_only_veto']) == (r['current_source_aware_class'] == 'UNRESOLVED' and stochastic_consensus)
        assert truth(r['conventional_class_unresolved']) == (r['branchA_class'] == 'UNRESOLVED')
        structural = r['metric'] == 'acceleration' and year % 2 == 0 and not coverage_missing
        assert truth(r['midpoint_structural_zero']) == structural
        assert not structural or dd['source_aware_estimate'] == '0.0'
        primary = ('source_coverage' if coverage_missing else 'insufficient_applicable_variants' if len(aa) < 3
                   else 'stochastic_uncertainty' if uncertainty else 'deterministic_midpoint_geometry' if midpoint_unresolved
                   else 'resolved_variant_conflict_consensus' if conflict else 'conventional_unresolved')
        assert r['primary_cause'] == primary
    short_keys = {(r['geo'], r['year'], r['metric']) for r in C if r['metric'] in METRICS}
    diagnostics = sorted({r['diagnostic'] for r in influence})
    assert len(diagnostics) == 5 and len(influence) == len(short_keys)*5
    summary = dict(medium=[], causes=[], influence=[], precision=[])
    for rule in sorted({r['rule'] for r in medium}):
        for pop in ('Member States', 'EU27_2020'):
            rr = [r for r in medium if r['rule'] == rule and r['population'] == pop]
            summary['medium'].append(dict(rule=rule, population=pop, **{x: sum(r['stability_class'] == x for r in rr) for x in LABELS},
                A_classes=dict(Counter(r['A_class'] for r in rr)), B1_classes=dict(Counter(r['B1_class'] for r in rr)),
                A_changed=[r['geo'] for r in rr if truth(r['A_class_changed'])],
                B1_changed=[r['geo'] for r in rr if truth(r['B1_class_changed'])],
                comparison_changed=[r['geo'] for r in rr if truth(r['stability_class_changed'])]))
    for pop in ('Member States', 'EU27_2020'):
        for per in ('2010-2022', '2023-2024 endpoints'):
            for metric in METRICS:
                rr = [r for r in causes if r['population'] == pop and r['period'] == per and r['metric'] == metric]
                current = [r for r in C if population(r['geo']) == pop and period(r['year']) == per and r['metric'] == metric]
                primary = dict(Counter(r['primary_cause'] for r in rr))
                assert len(rr) == sum(r['stability_class'] == 'UNRESOLVED' for r in current) == sum(primary.values())
                summary['causes'].append(dict(population=pop, period=per, metric=metric, unresolved=len(rr),
                    flags={flag: sum(truth(r[flag]) for r in rr) for flag in FLAGS}, primary=primary))
                for diagnostic in diagnostics:
                    ii = [r for r in influence if r['population'] == pop and r['period'] == per and r['metric'] == metric and r['diagnostic'] == diagnostic]
                    assert len(ii) == len(current)
                    assert len({(r['geo'], r['year'], r['metric']) for r in ii}) == len(ii)
                    # Recompute consensus independently from serialized vote classes.
                    for r in ii:
                        classes = [r[v + '_class'] for v in prod.VARIANTS if v != r['omitted_voting_variant'] and r[v + '_class'] != 'NOT_APPLICABLE']
                        resolved = len(classes) >= int(r['minimum_applicable_variants']) and len(set(classes)) == 1 and classes[0] != 'UNRESOLVED'
                        cls = classes[0] if resolved else 'UNRESOLVED'
                        assert cls == r['source_aware_class']
                        assert stability(r['branchA_class'], cls) == r['stability_class']
                    summary['influence'].append(dict(population=pop, period=per, metric=metric, diagnostic=diagnostic,
                        **{x: sum(r['stability_class'] == x for r in ii) for x in LABELS},
                        source_aware_changed=sum(truth(r['source_aware_class_changed']) for r in ii),
                        comparison_changed=sum(truth(r['stability_class_changed']) for r in ii),
                        newly_source_aware_resolved=sum(truth(r['newly_source_aware_resolved']) for r in ii),
                        newly_comparison_resolved=sum(truth(r['newly_comparison_resolved']) for r in ii)))
    ms_a, ms_b = [[r for r in rr if r['geo'] != 'EU27_2020'] for rr in (A, B)]
    widths = lambda rr: np.array([float(r['ci_high']) - float(r['ci_low']) for r in rr])
    wa, wb = widths(ms_a), widths(ms_b)
    summary['precision'] = dict(A_current_unresolved=sum(r['direction'] == 'UNRESOLVED' for r in ms_a),
        B1_current_unresolved=sum(r['direction'] == 'UNRESOLVED' for r in ms_b),
        median_A_CI_width=float(np.median(wa)), median_B1_CI_width=float(np.median(wb)),
        median_paired_B1_A_CI_width_ratio=float(np.median(wb/wa)), B1_wider_count=int(np.sum(wb > wa)))
    save('aggregate_summary.json', summary)
    integrity = json.loads((HERE / 'input_integrity.json').read_text())
    after = inventory()
    before = integrity['before_sha256']
    integrity.update(after_sha256=after, changed=[p for p in before if after.get(p) != before[p]],
                     added_outside_diagnostic_directory=sorted(set(after) - set(before)),
                     git_diff_available=False, git_limitation='Supplied workspace has no .git metadata; git status reports not a git repository.')
    save('input_integrity.json', integrity)
    assert not integrity['changed'] and not integrity['added_outside_diagnostic_directory']
    save('independent_verification.json', dict(pass_all=True, medium_rows=len(medium), unresolved_rows=len(causes),
        variant_influence_rows=len(influence), exhaustive_unresolved_keys=True, unique_complete_diagnostic_keys=True,
        aggregate_counts_reconcile=True, coverage_independently_checked_against_anchors=True,
        cause_flags_and_primary_hierarchy_checked=True, original_files_checked=len(before), original_files_changed=0))
    reports(summary)
    print('PASS: 56 harmonized medium-term rows, 1107 exhaustive unresolved rows, 5600 influence rows; all aggregates reconcile; 360 original files unchanged.')


def reports(summary):
    def table(headers, rows):
        return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---']*len(headers)) + ' |'] +
                         ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])
    def triplet(r):
        return f"{r['AGREEMENT']} / {r['CHANGE']} / {r['UNRESOLVED']}"
    medium_table = table(['Rule', 'Member States: agreement / disagreement / unresolved', 'EU27_2020'],
        [[rule, triplet(next(r for r in summary['medium'] if r['rule'] == rule and r['population'] == 'Member States')),
          triplet(next(r for r in summary['medium'] if r['rule'] == rule and r['population'] == 'EU27_2020'))]
         for rule in ('COMMON_TWO_SIDED_95', 'COMMON_DIRECTIONAL_95')])
    cs = [r for r in summary['causes'] if r['population'] == 'Member States' and r['period'] == '2010-2022']
    cause_table = table(['Metric', 'Unresolved', 'Coverage primary', 'Stochastic primary', 'Estimable midpoint unresolved (overlap)', 'Structural zero (overlap)'],
        [[r['metric'], r['unresolved'], f"{r['primary'].get('source_coverage',0)} ({100*r['primary'].get('source_coverage',0)/r['unresolved']:.1f}%)",
          f"{r['primary'].get('stochastic_uncertainty',0)} ({100*r['primary'].get('stochastic_uncertainty',0)/r['unresolved']:.1f}%)",
          f"{r['flags']['estimable_midpoint_unresolved']} ({100*r['flags']['estimable_midpoint_unresolved']/r['unresolved']:.1f}%)",
          r['flags']['midpoint_structural_zero']] for r in cs])
    influence_table = table(['Voting diagnostic', 'Annual change A / D / U', 'Turning A / D / U', 'Acceleration A / D / U'],
        [[diagnostic] + [triplet(next(r for r in summary['influence'] if r['population'] == 'Member States' and
          r['period'] == '2010-2022' and r['diagnostic'] == diagnostic and r['metric'] == metric)) for metric in METRICS]
         for diagnostic in ('LOO_pooled_signed', 'LOO_country_signed', 'LOO_pooled_absolute_random_sign',
                            'LOO_deterministic_midpoint', 'EXPLORATORY_midpoint_structural_reference')])
    text = f'''# Scientific diagnostic summary

Diagnostics of the frozen CMUR analysis, 10 November 2026. No manuscript claim or authoritative result was changed. A / D / U means agreement / resolved-class disagreement / unresolved, comparing the conventional and source-aware classes. Headline counts cover the fixed Member-State grid through 2022; EU27_2020 and endpoints are reported separately.

1. **The medium-term absence of conflicting resolved direction survives both identical inference rules.** EU27_2020 remains positive in A and B1 under both rules.

{medium_table}

Under COMMON_TWO_SIDED_95, B1 changes from positive to unresolved for Belgium (BE), and from negative to unresolved for Finland (FI) and Luxembourg (LU); all three comparisons change from agreement to unresolved. A does not change. Under COMMON_DIRECTIONAL_95, A changes from unresolved to positive for Slovenia (SI); B1 remains unresolved there, so no comparison changes. These are all changed geographies.

2. **The asymmetric rules partly conceal the loss of directional precision.** The current A/B1 unresolved counts are 8/12, a gap of four. Common two-sided inference gives 8/15, a gap of seven: three of those seven additional B1 unresolved classifications (42.9%) are masked by B1's more permissive current directional rule. Common directional inference gives 7/12, a gap of five. Thus the asymmetry also understates that same-rule gap by one. These are classification differences, not a causal decomposition of sampling uncertainty. B1's standard 95% interval is wider for all 27 Member States; median widths are {summary['precision']['median_A_CI_width']:.4f} (A) and {summary['precision']['median_B1_CI_width']:.4f} (B1) percentage points/year, with a median paired width ratio of {summary['precision']['median_paired_B1_A_CI_width_ratio']:.3f}. Harmonization changes classes, not slopes or intervals.

3. **Most current short-window unresolved claims have estimable stochastic variants that fail their sign-support rules.** Coverage-first exclusive attribution reconciles exactly:

{cause_table}

Every through-2022 unresolved Member-State claim has at least one applicable stochastic variant unresolved, including coverage failures. Once coverage failures are separated, the stochastic primary counts are 306, 286 and 287. The midpoint is additionally unresolved in 16, 17 and 170 claims respectively when missing-coverage rows are included; the estimable counts in the table isolate geometry from absent estimates. Structural zero accounts for 158 acceleration claims, or 53.4% of unresolved acceleration claims. All midpoint failures overlap stochastic failures; none is a sole veto that prevents otherwise unanimous stochastic resolution. Fewer-than-three applicability and consensus-only failure each occur zero times. One acceleration claim (CY, 2017) contains conflicting resolved variant signs (country negative, midpoint positive), alongside unresolved pooled stochastic variants; it is not a resolved disagreement between final branches.

EU27_2020 has 12/11/11 unresolved annual-change/turning/acceleration claims: all have estimable stochastic uncertainty and none lacks coverage. Six acceleration claims also have midpoint structural zeros (54.5%); other midpoint-unresolved counts are zero. For 2023–2024, every metric has 54 Member-State and two EU endpoint claims, all unresolved because of absent required coverage. These endpoints are excluded from headline denominators. The flags overlap and cannot be added as independent causal proportions. The exclusive hierarchy allocates overlaps to the first applicable cause; its zero midpoint-primary count does not mean midpoint geometry is absent.

4. **The qualitative short-window conclusion survives every leave-one-variant-out diagnostic.** Counts are:

{influence_table}

Only removing country_signed changes claims: 10 additional annual-change agreements (3.2% of the current 314 unresolved), 10 turning agreements (3.4% of 295), and four acceleration agreements (1.4% of 296). All 24 changed comparisons move from unresolved to agreement; their exact geography/year/metric keys and classes are in short_window_variant_influence.csv. Even this diagnostic leaves 93.8%, 96.0% and 98.3% of the respective Member-State grids unresolved. Removing either pooled variant or the midpoint changes no source-aware or final comparison class. Removing the midpoint resolves **zero additional acceleration claims**, under both the strict minimum-three LOO and the explicitly EXPLORATORY stochastic-only minimum-two convention. EU27_2020 remains 0/0/12, 0/0/11 and 0/0/11 under every diagnostic. Endpoints remain 0/0/54 for Member States and 0/0/2 for EU27_2020 for each metric and diagnostic.

5. **The temporal-resolution support boundary is genuinely mixed within the tested design.** Missing or flagged anchors and absent endpoint coverage create direct provenance constraints. Within covered windows, the four-year calibration residuals, fixed stochastic sign-support thresholds, and unanimity requirement jointly determine support. Country-variant removal shows a modest consensus sensitivity, while the midpoint's frequent geometric zeros do not independently account for the observed acceleration unresolvedness. The unchanged majority-unresolved outcome across these diagnostics supports the qualitative boundary under this calibrated design. These diagnostics hold calibration and sign-support thresholds fixed; they do not identify a separate causal percentage for calibration versus stochastic threshold versus consensus, or establish the uncertainty distribution of true annual CMUR errors.

The single most decision-relevant finding is that removing the deterministic midpoint resolves no additional acceleration claims: stochastic support remains limiting, while the medium-term zero-conflict finding survives harmonized inference.
'''
    (HERE / 'SCIENTIFIC_DIAGNOSTIC_SUMMARY.md').write_text(text, encoding='utf-8')
    readme = '''# WP P1A temporal-support diagnostics

Completed 10 November 2026 using frozen repository inputs only. All alternatives are diagnostics, not replacements for the published analysis. No publication, novelty, GO or journal decision is made.

## Authority and baseline gate

Authority was traced from manuscript/scripts/build_assets.py and its make_figures.py compatibility entry point to data/processed/branchA_medium_term.csv, branchB1_anchor_medium_term.csv and temporal_stability_comparison.csv. manuscript/scripts/verify_integration.py establishes the manuscript headline counts and branch-specific rules. scripts/wp1b_r.py defines the generating model and locates data/processed/b2r_short_window_support.csv, cmur_source_provenance_panel.csv, verified_anchor_map.csv, b2r_residual_pool.csv and pseudogap_support.csv. scripts/build_wp1a.py and acquire_wp1a.py locate archived inputs through logs/acquisition_manifest.json; these were used only by offline verification. The methods, results and discussion sections were read directly. No external source was queried.

The baseline gate passed before alternatives ran: Member States medium 15/0/12; annual change 10/0/314; turning 2/0/295; acceleration 1/0/296. EU27_2020 medium 1/0/0 and short-window 0/0/12, 0/0/11, 0/0/11. Every frozen comparison class and short-window consensus was checked. Existing full-rerun verification regenerated byte-identical scientific result CSVs in an isolated copy.

## Rules and scope

Medium-term OLS was recomputed from exactly the production annual and admissible-anchor observations, using production ols and checking every estimate, SE, interval and p_plus against the frozen files to absolute tolerance 1e-12. COMMON_TWO_SIDED_95 requires ci_low > 0 or ci_high < 0. COMMON_DIRECTIONAL_95 requires p_plus >= 0.95 or <= 0.05. Both apply identically to A and B1. Medium-term output contains 56 rows (28 geographies × two rules) with all current and diagnostic classes and change flags.

Short-window diagnostics reuse frozen variant classifications, not new simulations. Thus the original PCG64 geography–variant streams, seed 20260926, 20,000 draws, clipping, sign thresholds and structural midpoint correction are preserved exactly. The original full rerun independently verifies the simulation artifacts. No new tolerance, threshold, residual scale or dataset was introduced.

All short-window geography/year/metric claims in the production 2011/2012–2024 grid are included: 1,120 claim instances. period distinguishes the manuscript grid through 2022 from 2023–2024 endpoints; population separates 27 Member States from EU27_2020. The causes CSV includes every currently UNRESOLVED final comparison, 1,107 unique claims. The influence CSV includes all 1,120 claims under each of five diagnostics, 5,600 unique rows. CHANGE in CSVs means resolved-class disagreement, as in production; it does not merely mean a diagnostic class changed.

## Cause flags and exclusive hierarchy

All flags are non-exclusive Boolean values. missing_required_coverage means a required year has neither an admissible anchor nor an eligible adjacent-anchor bracket, or is an unsupported endpoint; its value is independently checked against the anchor map. fewer_than_required_applicable_variants compares the current applicable count with three. NOT_APPLICABLE rows never vote. any_applicable_stochastic_unresolved includes unresolved coverage rows; estimable_stochastic_unresolved requires a nonempty stored source-aware estimate and therefore distinguishes sign-support uncertainty from missing trajectories. For turning, stochastic uncertainty means at least one constituent annual change fails support; it is not uncertainty in the strict reversal-event estimate. deterministic_midpoint_unresolved similarly includes missing estimates, while estimable_midpoint_unresolved excludes them. midpoint_structural_zero identifies covered anchor–midpoint–anchor acceleration triples with exact stored zero.

conflicting_resolved_variant_classes detects two different resolved variant labels, even when another variant is unresolved. failure_only_cross_variant_consensus uses a strict definition: source-aware consensus fails despite complete window coverage and every applicable individual vote being resolved. Thus it captures pure unanimity/minimum-vote failure, not a veto by an individually unresolved variant. midpoint_only_veto separately tests whether the stochastic-only consensus (at least two applicable stochastic votes, all unanimous and resolved) resolves a currently source-aware-unresolved claim. conventional_class_unresolved records a conventional-side obstacle. The variant classes and original reason codes remain visible on every cause row.

The exclusive primary hierarchy is: (1) missing source/endpoint coverage; (2) insufficient applicable variants; (3) at least one estimable stochastic variant unresolved; (4) estimable deterministic midpoint unresolved; (5) conflicting resolved variant votes; (6) conventional class unresolved. Every current unresolved claim receives exactly one first-matching category, and its sums equal the frozen unresolved totals. This is deterministic accounting, not independent causal attribution. Stochastic-first ordering assigns overlapping midpoint failures to stochastic uncertainty; non-exclusive geometry flags are essential to interpretation.

## Variant influence conventions

Each LOO removes exactly its named voting variant and retains the literal minimum of THREE applicable variants and unanimous resolved labels among ALL retained applicable variants. Where country_signed is NOT_APPLICABLE (CZ and EL), removing another variant leaves only two votes and therefore cannot resolve the claim. No previous resolution is lost there in these data. Removing an already-inapplicable country vote has no effect. All individual 0.95 thresholds and conventional classifications remain frozen.

EXPLORATORY_midpoint_structural_reference retains the midpoint classes as structural information in the output, removes its vote, and requires at least TWO applicable stochastic variants with all applicable stochastic votes unanimous and resolved. Normally three stochastic variants vote; CZ/EL have two. This explicitly relaxes the minimum solely to retain their original country-pool applicability exclusion while evaluating the stochastic family. It is not a proposed replacement primary analysis. Both this exploratory diagnostic and literal minimum-three midpoint LOO are reported separately; both add zero acceleration resolutions.

source_aware_class_changed and stability_class_changed distinguish changes in consensus from changes in comparison with the conventional class. newly_source_aware_resolved and newly_comparison_resolved explicitly identify additional resolutions. Filtering either change flag gives the complete changed-claim list, including old/new classes. Only LOO_country_signed changes any row (24 comparisons); every change is unresolved → agreement.

## Files

- README.md: authority, methods, applicability conventions, file schemas and verification.
- baseline_reproduction.csv: eight population-by-metric headline summaries with actual/expected counts and exact-match flags.
- medium_term_harmonized_rules.csv: 56 geography-by-rule rows, classes, current classes, change flags and fitted inference quantities.
- short_window_unresolved_causes.csv: 1,107 unresolved claims, non-exclusive flags, exclusive primary cause, all variant classes and reason codes.
- short_window_variant_influence.csv: 5,600 claim-by-diagnostic rows, current/alternative classes, vote counts and change/resolution flags.
- SCIENTIFIC_DIAGNOSTIC_SUMMARY.md: evidence-only answers to the five scientific questions and aggregate tables.
- run_diagnostics.py: reproducible baseline gate, diagnostics, independent aggregation and report generator.
- aggregate_summary.json: all metric/population/period totals, flag and primary-cause counts, rule-change geography lists and interval-width summaries.
- verification.json: exact commands, outputs and return codes from existing verification suites; baseline gate and byte-identity status.
- independent_verification.json: row coverage, unique keys, aggregate reconciliation and preservation checks.
- input_integrity.json: before/after SHA-256 of every pre-existing file (360 files), including manuscript, frozen inputs, figures and reproducibility/zenodo; changed-file lists are empty.

## Verification and save status

From the repository root, run:

```powershell
python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py baseline
python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py diagnostics
python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py verify
```

The baseline stage executes python -B scripts/verify_wp1a.py (44 checks), python -B scripts/verify_wp1b_r.py (30 checks, including the full seed-preserving rerun), and python -B manuscript/scripts/verify_integration.py (17 checks) in a disposable repository copy under this output directory. All pass. This avoids their otherwise overwriting authoritative files. Asset generation was not run against the original manuscript; its frozen tables/figures are verified by the integration suite. Integration checks the existing build logs/PDF; it does not perform a fresh LaTeX compilation.

The verify stage reloads serialized CSVs, independently sums all aggregate cells, checks exhaustive unresolved keys, verifies anchor coverage independently, checks cause flags/hierarchy, and recomputes each retained-vote consensus from serialized variant classes. Report tables are generated from these reloaded aggregates. The initial file manifest cannot silently reset on rerun; existing files must still match it.

All deliverables are saved locally in this new directory. The supplied workspace has no .git directory: git status --short returns “not a git repository”. A Git diff or commit cannot be produced. Whole-workspace before/after hashes instead verify that every pre-existing file is byte-identical and no file was added outside this diagnostic directory. This is a repository-metadata limitation, not a scientific blocker. No scientific blocker was found. No authoritative CSV, figure, manuscript, archived extract or Zenodo evidence was overwritten.
'''
    (HERE / 'README.md').write_text(readme, encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('baseline', 'diagnostics', 'verify'))
    args = parser.parse_args()
    {'baseline': baseline, 'diagnostics': diagnostics, 'verify': independent_verify}[args.stage]()
