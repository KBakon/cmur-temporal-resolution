"""Paired amplitude-only CMUR sensitivity diagnostics. No writes outside HERE."""
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
from collections import Counter, defaultdict
from decimal import Decimal

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / 'analysis/wp_p1a_temporal_support_diagnostics'
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'scripts'))
import wp1b_r as P

U = 'UNSCALED_STRESS'
B = 'BROWNIAN_SD_1_SQRT2'
SCALES = {U: 1.0, B: float(1 / np.sqrt(2))}
GRID = (.80, .85, .90, .925, .95, .975, .99)
METRICS = ('annual_change', 'turning', 'acceleration')
POPS = ('Member States', 'EU27_2020')
STOCH = P.VARIANTS[:-1]


def read(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def write(name, rows, fields=None):
    with (HERE / name).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), lineterminator='\n')
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


def array_hash(x):
    return hashlib.sha256(np.asarray(x).tobytes()).hexdigest()


def inventory():
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(ROOT.rglob('*'))
            if p.is_file() and not p.is_relative_to(HERE)}


def pop(geo):
    return POPS[1] if geo == 'EU27_2020' else POPS[0]


def compare(a, b):
    return 'UNRESOLVED' if 'UNRESOLVED' in (a, b) else 'AGREEMENT' if a == b else 'CHANGE'


def baseline():
    initial = inventory()
    if (HERE / 'input_integrity.json').exists():
        assert initial == json.loads((HERE / 'input_integrity.json').read_text())['before_sha256']
    else:
        save('input_integrity.json', dict(before_sha256=initial))
    # All mandatory prior evidence is read and checked before any alternative.
    prior = {name: read(PRIOR / name) for name in ('baseline_reproduction.csv', 'medium_term_harmonized_rules.csv',
             'short_window_unresolved_causes.csv', 'short_window_variant_influence.csv')}
    prior_json = {name: json.loads((PRIOR / name).read_text()) for name in
                  ('aggregate_summary.json', 'verification.json', 'independent_verification.json')}
    assert (PRIOR / 'SCIENTIFIC_DIAGNOSTIC_SUMMARY.md').read_text(encoding='utf-8')
    assert prior_json['verification.json']['baseline_pass'] and prior_json['independent_verification.json']['pass_all']
    C = read(P.OUT / 'temporal_stability_comparison.csv')
    gate = []
    expected = {'medium_term_slope': (15, 0, 12), 'annual_change': (10, 0, 314),
                'turning': (2, 0, 295), 'acceleration': (1, 0, 296)}
    for metric, target in expected.items():
        for population in POPS:
            rr = [r for r in C if r['metric'] == metric and pop(r['geo']) == population and
                  (metric == 'medium_term_slope' or int(r['year']) <= 2022)]
            actual = tuple(sum(r['stability_class'] == cls for r in rr) for cls in ('AGREEMENT', 'CHANGE', 'UNRESOLVED'))
            wanted = target if population == POPS[0] else (1, 0, 0) if metric == 'medium_term_slope' else (0, 0, 12 if metric == 'annual_change' else 11)
            assert actual == wanted
            old = next(r for r in prior['baseline_reproduction.csv'] if r['metric'] == metric and r['population'] == population)
            assert actual == tuple(int(old[k]) for k in ('agreement', 'resolved_class_disagreement', 'unresolved'))
            gate.append(dict(check='production_' + metric, population=population, observed='/'.join(map(str, actual)), expected='/'.join(map(str, wanted)), passed=True))
    frozen_a = {r['geo']: r for r in read(P.OUT / 'branchA_medium_term.csv')}
    frozen_b = {r['geo']: r for r in read(P.OUT / 'branchB1_anchor_medium_term.csv')}
    for rule, target in (('COMMON_TWO_SIDED_95', (12, 0, 15)), ('COMMON_DIRECTIONAL_95', (15, 0, 12))):
        for population in POPS:
            rr = [r for r in prior['medium_term_harmonized_rules.csv'] if r['rule'] == rule and r['population'] == population]
            for r in rr:
                cls = []
                for branch, frozen in (('A', frozen_a[r['geo']]), ('B1', frozen_b[r['geo']])):
                    label = ('POSITIVE' if float(frozen['ci_low']) > 0 else 'NEGATIVE' if float(frozen['ci_high']) < 0 else 'UNRESOLVED') if rule == 'COMMON_TWO_SIDED_95' else ('POSITIVE' if float(frozen['p_positive']) >= .95 else 'NEGATIVE' if float(frozen['p_positive']) <= .05 else 'UNRESOLVED')
                    assert label == r[branch + '_class']
                    cls.append(label)
                assert compare(*cls) == r['stability_class']
            actual = tuple(sum(r['stability_class'] == cls for r in rr) for cls in ('AGREEMENT', 'CHANGE', 'UNRESOLVED'))
            wanted = target if population == POPS[0] else (1, 0, 0)
            assert actual == wanted
            reported = next(r for r in prior_json['aggregate_summary.json']['medium'] if r['rule'] == rule and r['population'] == population)
            assert actual == tuple(reported[k] for k in ('AGREEMENT', 'CHANGE', 'UNRESOLVED'))
            gate.append(dict(check=rule, population=population, observed='/'.join(map(str, actual)), expected='/'.join(map(str, wanted)), passed=True))
    frozen_s = read(P.OUT / 'b2r_short_window_support.csv')
    lookup = defaultdict(list)
    for r in frozen_s:
        if r['variant'] != 'deterministic_midpoint' and r['variant_class'] != 'NOT_APPLICABLE':
            lookup[(r['geo'], r['year'], r['metric'])].append(r['variant_class'])
    for diagnostic, minimum in (('LOO_deterministic_midpoint', 3), ('EXPLORATORY_midpoint_structural_reference', 2)):
        for population in POPS:
            rr = [r for r in prior['short_window_variant_influence.csv'] if r['diagnostic'] == diagnostic and
                  r['metric'] == 'acceleration' and r['population'] == population and int(r['year']) <= 2022]
            for r in rr:
                cc = lookup[(r['geo'], r['year'], r['metric'])]
                final = cc[0] if len(cc) >= minimum and len(set(cc)) == 1 and cc[0] != 'UNRESOLVED' else 'UNRESOLVED'
                assert final == r['source_aware_class']
            newly = sum(r['newly_source_aware_resolved'] == 'True' for r in rr)
            assert newly == 0
            gate.append(dict(check=diagnostic + '_new_acceleration_resolutions', population=population, observed=str(newly), expected='0', passed=True))
    cause_keys = {(r['geo'], r['year'], r['metric']) for r in prior['short_window_unresolved_causes.csv']}
    assert cause_keys == {(r['geo'], r['year'], r['metric']) for r in C if r['metric'] in METRICS and r['stability_class'] == 'UNRESOLVED'}
    verify_prior_details()
    write('baseline_gate.csv', gate)
    print('PASS production and P1A gates; starting isolated production verification.', flush=True)
    sandbox = (HERE / 'verification_workspace').resolve()
    assert sandbox.is_relative_to(HERE) and sandbox != HERE
    sandbox.mkdir(exist_ok=False)  # avoid Windows tempfile mode-0700 ACL issue
    suites = []
    try:
        for folder in ('scripts', 'data', 'logs', 'manuscript'):
            shutil.copytree(ROOT / folder, sandbox / folder, ignore=shutil.ignore_patterns('__pycache__'))
        for script in ('scripts/verify_wp1a.py', 'scripts/verify_wp1b_r.py', 'manuscript/scripts/verify_integration.py'):
            result = subprocess.run([sys.executable, '-B', script], cwd=sandbox, capture_output=True, text=True,
                                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
            suites.append(dict(command='python -B ' + script, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
            save('verification.json', dict(baseline_pass=False, existing_suites=suites))
            assert result.returncode == 0, suites[-1]
        outputs = __import__('verify_wp1a').PRODUCTS + __import__('verify_wp1b_r').FILES
        assert all(sha(P.OUT / n) == sha(sandbox / 'data/processed' / n) for n in outputs)
    finally:
        assert sandbox.is_relative_to(HERE) and sandbox != HERE
        shutil.rmtree(sandbox)
    assert inventory() == initial
    save('verification.json', dict(baseline_pass=True, existing_suites=suites,
        baseline_regenerated_csv_bytes_identical=True, seed=P.SEED, trajectories=P.N_DRAWS))
    print('PASS all existing suites and byte-identical frozen production rerun. Alternatives authorized by gate.', flush=True)


def paired_paths(panel, anchors, by_geo, pooled, eligible, geo, variant):
    if variant == 'country_signed' and len(by_geo[geo]) < 4:
        return {U: None, B: None}, []
    det = variant == 'deterministic_midpoint'
    n = 1 if det else P.N_DRAWS
    gi, vi = P.GEOS.index(geo), P.VARIANTS.index(variant)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([P.SEED, gi, vi])))
    # One identity/sign sequence feeds both scenarios; independent replay below checks it.
    replay = np.random.Generator(np.random.PCG64(np.random.SeedSequence([P.SEED, gi, vi])))
    paths = {U: {}, B: {}}
    audit = []
    for year in P.YEARS:
        if P.admissible(anchors, geo, year):
            fixed = np.full(n, P.cmur(panel, geo, year), dtype=float)
            paths[U][year] = fixed.copy()
            paths[B][year] = fixed.copy()
        elif (geo, year) in eligible:
            midpoint = eligible[(geo, year)]
            if det:
                raw = np.zeros(n)
                ids = np.zeros(n, dtype=np.int64)
                signs = np.ones(n)
            else:
                pool = np.asarray(by_geo[geo], dtype=float) if variant == 'country_signed' else pooled
                if variant == 'pooled_absolute_random_sign':
                    pool = np.abs(pool)
                ids = rng.choice(len(pool), size=n, replace=True)
                signs = rng.choice(np.array([-1.0, 1.0]), size=n, replace=True) if variant == 'pooled_absolute_random_sign' else np.ones(n)
                raw = pool[ids] * signs
                # Replay production's choice-of-values API; array equality also protects draw order.
                check = replay.choice(pool, size=n, replace=True)
                if variant == 'pooled_absolute_random_sign':
                    check *= replay.choice(np.array([-1.0, 1.0]), size=n, replace=True)
                assert np.array_equal(raw, check)
            scaled = raw * SCALES[B]
            assert np.array_equal(scaled, raw * SCALES[B])
            paths[U][year] = np.clip(midpoint + raw, 0, 100)
            paths[B][year] = np.clip(midpoint + scaled, 0, 100)
            audit.append(dict(geo=geo, variant=variant, year=year, draws=n, stream_seed=P.SEED,
                geography_stream_index=gi, variant_stream_index=vi, midpoint=midpoint,
                U_identity_sha256=array_hash(ids), B_identity_sha256=array_hash(ids),
                U_sign_sha256=array_hash(signs), B_sign_sha256=array_hash(signs),
                U_residual_sha256=array_hash(raw), B_residual_sha256=array_hash(scaled),
                B_scale=SCALES[B], exact_amplitude_transform=True, independent_stream_replay_equal=True,
                U_clipped_draws=int(np.sum((midpoint+raw < 0) | (midpoint+raw > 100))),
                B_clipped_draws=int(np.sum((midpoint+scaled < 0) | (midpoint+scaled > 100)))))
    original = P.simulate(panel, anchors, by_geo, pooled, eligible, geo, variant)
    assert set(paths[U]) == set(original)
    assert all(np.array_equal(paths[U][y], original[y]) for y in original), (geo, variant)
    return paths, audit


def variant_rows(panel, eligible, geo, variant, path):
    """Reproduce every original field, including production row order and turning reasons."""
    rows, annual = [], {}
    det = variant == 'deterministic_midpoint'
    def calculate(year, metric):
        required = (year-1, year) if metric == 'annual_change' else (year-2, year-1, year)
        defined = path is not None and year <= 2022 and all(y in path for y in required)
        if not defined:
            cls = 'NOT_APPLICABLE' if path is None else 'UNRESOLVED'
            reason = 'country_pool_lt_4' if path is None else 'endpoint_unanchored' if year >= 2023 else 'missing_admissible_anchor_or_bracket'
            result = dict(estimate='', low='', high='', p_positive='', p_negative='', classification=cls, sign_support='')
        elif metric != 'turning':
            values = np.zeros(1, dtype=float) if metric == 'acceleration' and det and (geo, year-1) in eligible else P.short_metric(path, year, metric)
            result = P.classify_array(values, det)
            reason = ''
        else:
            a1, a2 = annual[year-1], annual[year]
            c1, c2 = a1['variant_class'], a2['variant_class']
            cls = 'UNRESOLVED' if 'UNRESOLVED' in (c1, c2) else 'TURN' if c1 != c2 else 'NO_TURN'
            result = dict(estimate=float(np.mean(P.short_metric(path, year, metric))), low='', high='',
                p_positive='', p_negative='', classification=cls,
                sign_support=min(a1['sign_support'], a2['sign_support']) if not det else '')
            reason = '' if cls != 'UNRESOLVED' else 'consecutive_change_sign_support_below_0.95'
        return dict(geo=geo, year=year, metric=metric, variant=variant,
            branchA_estimate=P.conventional_estimate(panel, geo, year, metric), branchA_class=P.conventional_class(panel, geo, year, metric),
            source_aware_estimate=result['estimate'], ci_low=result['low'], ci_high=result['high'],
            p_positive=result['p_positive'], p_negative=result['p_negative'], sign_support=result['sign_support'],
            variant_class=result['classification'], reason=reason)
    for year in range(2011, 2025):
        row = calculate(year, 'annual_change')
        rows.append(row)
        annual[year] = row
    for year in range(2012, 2025):
        for metric in ('turning', 'acceleration'):
            rows.append(calculate(year, metric))
    return rows, annual


def preferred(pp, pn):
    return 'POSITIVE' if pp > pn else 'NEGATIVE' if pn > pp else 'TIE'


def continuous(row, annual, scenario):
    result = dict(scenario=scenario, residual_amplitude_scale=SCALES[scenario], population=pop(row['geo'])) | row.copy()
    extra = dict(applicable=row['variant_class'] != 'NOT_APPLICABLE', estimable=row['source_aware_estimate'] != '',
        stochastic=row['variant'] != 'deterministic_midpoint', preferred_sign='', preferred_class='',
        preferred_sign_support='', distance_from_0_95_threshold='', support_definition='',
        earlier_change_p_positive='', earlier_change_p_negative='', earlier_change_preferred_sign='', earlier_change_support='', earlier_change_class_0_95='',
        later_change_p_positive='', later_change_p_negative='', later_change_preferred_sign='', later_change_support='', later_change_class_0_95='')
    if extra['stochastic'] and extra['estimable']:
        if row['metric'] != 'turning':
            pp, pn = row['p_positive'], row['p_negative']
            extra.update(preferred_sign=preferred(pp, pn), preferred_class=preferred(pp, pn),
                         preferred_sign_support=max(pp, pn), support_definition='max(p_positive,p_negative)')
        else:
            prefs, supports = [], []
            for prefix, year in (('earlier', row['year']-1), ('later', row['year'])):
                aa = annual[year]
                pp, pn = aa['p_positive'], aa['p_negative']
                pref, support = preferred(pp, pn), max(pp, pn)
                extra.update({prefix+'_change_p_positive': pp, prefix+'_change_p_negative': pn,
                              prefix+'_change_preferred_sign': pref, prefix+'_change_support': support,
                              prefix+'_change_class_0_95': aa['variant_class']})
                prefs.append(pref)
                supports.append(support)
            cls = 'TIE' if 'TIE' in prefs else 'TURN' if prefs[0] != prefs[1] else 'NO_TURN'
            extra.update(preferred_class=cls, preferred_sign_support=min(supports),
                         support_definition='minimum_of_two_constituent_max_sign_supports_not_turn_event_probability')
        extra['distance_from_0_95_threshold'] = extra['preferred_sign_support'] - .95
    elif not extra['stochastic']:
        extra['support_definition'] = 'structural_reference_no_probability'
    return result | extra


def run_scenarios():
    verification = json.loads((HERE / 'verification.json').read_text())
    assert verification['baseline_pass']
    assert inventory() == json.loads((HERE / 'input_integrity.json').read_text())['before_sha256']
    panel, anchors = P.load_panel()
    residuals, by_geo, pooled, eligible = P.prepare(panel, anchors)
    assert len(pooled) == 130
    frozen_r = read(P.OUT / 'b2r_residual_pool.csv')
    assert [str(r['signed_residual_pp']) for r in residuals] == [r['signed_residual_pp'] for r in frozen_r]
    all_rows, audit, u_original = [], [], []
    # Complete all U checks before constructing any B output/classification.
    paths_by = {}
    for geo in P.GEOS:
        for variant in P.VARIANTS:
            paths, aa = paired_paths(panel, anchors, by_geo, pooled, eligible, geo, variant)
            paths_by[(geo, variant)] = paths
            audit.extend(aa)
            rr, annual = variant_rows(panel, eligible, geo, variant, paths[U])
            u_original.extend(rr)
            all_rows.extend(continuous(r, annual, U) for r in rr if r['year'] <= 2022)
    frozen = read(P.OUT / 'b2r_short_window_support.csv')
    write('scenario_u_production_reproduction.csv', u_original, list(frozen[0]))
    assert sha(HERE / 'scenario_u_production_reproduction.csv') == sha(P.OUT / 'b2r_short_window_support.csv'), 'STOP: U is not byte-identical to production.'
    verification['scenario_U_production_csv_byte_identical'] = True
    save('verification.json', verification)
    print('PASS Scenario U production output byte identity; evaluating amplitude-only B.', flush=True)
    for geo in P.GEOS:
        for variant in P.VARIANTS:
            rr, annual = variant_rows(panel, eligible, geo, variant, paths_by[(geo, variant)][B])
            all_rows.extend(continuous(r, annual, B) for r in rr if r['year'] <= 2022)
    write('calibration_scenario_support.csv', all_rows)
    write('paired_draw_audit.csv', audit)
    verification.update(paired_identity_sign_order_verified=True, independent_stream_replay_equal=True,
                        U_paths_exactly_equal_production=True, only_residual_amplitude_changed=True,
                        brownian_amplitude_scale=SCALES[B], paired_draw_records=len(audit))
    save('verification.json', verification)
    del paths_by
    multispan(panel, anchors)
    aggregate()


def stochastic_label(r, threshold):
    if r['variant_class'] == 'NOT_APPLICABLE':
        return 'NOT_APPLICABLE'
    if r['source_aware_estimate'] == '':
        return 'UNRESOLVED'
    if r['variant'] == 'deterministic_midpoint':
        return r['variant_class']
    if r['metric'] != 'turning':
        return 'POSITIVE' if float(r['p_positive']) >= threshold else 'NEGATIVE' if float(r['p_negative']) >= threshold else 'UNRESOLVED'
    c = []
    for prefix in ('earlier', 'later'):
        c.append('POSITIVE' if float(r[prefix+'_change_p_positive']) >= threshold else
                 'NEGATIVE' if float(r[prefix+'_change_p_negative']) >= threshold else 'UNRESOLVED')
    return 'UNRESOLVED' if 'UNRESOLVED' in c else 'TURN' if c[0] != c[1] else 'NO_TURN'


def aggregate():
    """All aggregates/prose are derived from reloaded, serialized outputs."""
    rows = read(HERE / 'calibration_scenario_support.csv')
    assert len(rows) == 2*28*34*4
    grouped = defaultdict(dict)
    for r in rows:
        key = (r['scenario'], r['geo'], r['year'], r['metric'])
        assert r['variant'] not in grouped[key]
        grouped[key][r['variant']] = r
        if r['variant'] == 'deterministic_midpoint':
            assert not r['p_positive'] and not r['p_negative'] and not r['preferred_sign_support']
        elif r['source_aware_estimate'] != '':
            if r['metric'] == 'turning':
                assert not r['p_positive'] and not r['p_negative']
                support = min(float(r[x+'_change_support']) for x in ('earlier', 'later'))
            else:
                support = max(float(r['p_positive']), float(r['p_negative']))
            assert support == float(r['preferred_sign_support'])
            assert support - .95 == float(r['distance_from_0_95_threshold'])
        assert stochastic_label(r, .95) == r['variant_class']
    claims, thresholds = [], []
    for key, vv in grouped.items():
        scenario, geo, year, metric = key
        ss = [r for v, r in vv.items() if v in STOCH and r['variant_class'] != 'NOT_APPLICABLE']
        available = all(r['source_aware_estimate'] != '' for r in ss)
        prefs = {r['preferred_class'] for r in ss}
        status = 'MISSING_COVERAGE' if not available else 'PREFERRED_TIE' if 'TIE' in prefs else 'CONCORDANT' if len(prefs) == 1 else 'CONFLICT'
        weakest = min(float(r['preferred_sign_support']) for r in ss) if status == 'CONCORDANT' else ''
        base = dict(scenario=scenario, population=pop(geo), geo=geo, year=year, metric=metric)
        claim = base | dict(applicable_stochastic_variants=';'.join(r['variant'] for r in ss), n_applicable_stochastic=len(ss),
            preferred_class_status=status, preferred_class=next(iter(prefs)) if status == 'CONCORDANT' else '',
            weakest_variant_support=weakest, midpoint_structural_class=vv['deterministic_midpoint']['variant_class'],
            midpoint_estimate=vv['deterministic_midpoint']['source_aware_estimate'],
            branchA_class=vv['pooled_signed']['branchA_class'])
        for v in STOCH:
            claim[v+'_preferred_class'] = vv[v]['preferred_class']
            claim[v+'_preferred_support'] = vv[v]['preferred_sign_support']
        claims.append(claim)
        for threshold in GRID:
            votes = {v: stochastic_label(r, threshold) for v, r in vv.items()}
            stochastic_votes = [votes[r['variant']] for r in ss]
            resolved = set(stochastic_votes) - {'UNRESOLVED'}
            # Conflict precedence: incompatible resolved stochastic labels count as conflict even if a third is unresolved.
            st = 'CONFLICT' if len(resolved) > 1 else 'RESOLVED_CONCORDANT' if len(ss) >= 2 and 'UNRESOLVED' not in stochastic_votes else 'UNRESOLVED'
            applicable = [cls for cls in votes.values() if cls != 'NOT_APPLICABLE']
            final = applicable[0] if len(applicable) >= 3 and len(set(applicable)) == 1 and applicable[0] != 'UNRESOLVED' else 'UNRESOLVED'
            thresholds.append(base | dict(threshold=threshold, stochastic_only_status=st,
                stochastic_only_class=stochastic_votes[0] if st == 'RESOLVED_CONCORDANT' else '',
                source_aware_class=final, branchA_class=claim['branchA_class'], stability_class=compare(claim['branchA_class'], final)) |
                {v+'_class': votes[v] for v in P.VARIANTS})
    write('claim_level_continuous_support.csv', claims)
    write('threshold_claim_classifications.csv', thresholds)
    # Reload again before summing any counts intended for reporting.
    claims = read(HERE / 'claim_level_continuous_support.csv')
    thresholds = read(HERE / 'threshold_claim_classifications.csv')
    survival = []
    distributions = []
    common = {(r['geo'], r['year'], r['metric']) for r in claims if r['scenario'] == U and r['preferred_class_status'] == 'CONCORDANT'} & {
              (r['geo'], r['year'], r['metric']) for r in claims if r['scenario'] == B and r['preferred_class_status'] == 'CONCORDANT'}
    for scenario in SCALES:
        for population in POPS:
            for metric in METRICS:
                cc = [r for r in claims if r['scenario'] == scenario and r['population'] == population and r['metric'] == metric]
                for subset in ('scenario_concordant', 'common_concordant_pairs'):
                    sample = cc if subset == 'scenario_concordant' else [r for r in cc if (r['geo'], r['year'], r['metric']) in common]
                    vals = np.array([float(r['weakest_variant_support']) for r in sample if r['weakest_variant_support'] != ''])
                    q = np.quantile(vals, [.1, .25, .5, .75, .9]) if len(vals) else ['']*5
                    distributions.append(dict(scenario=scenario, population=population, metric=metric, subset=subset,
                        total_claims=len(cc), concordant_with_support=len(vals), conflict=sum(r['preferred_class_status'] == 'CONFLICT' for r in cc),
                        preferred_tie=sum(r['preferred_class_status'] == 'PREFERRED_TIE' for r in cc),
                        missing_coverage=sum(r['preferred_class_status'] == 'MISSING_COVERAGE' for r in cc),
                        p10=q[0], q25=q[1], median=q[2], q75=q[3], p90=q[4],
                        **{f'proportion_ge_{cut:.2f}': float(np.mean(vals >= cut)) if len(vals) else '' for cut in (.80, .90, .95, .99)}))
                for threshold in GRID:
                    tt = [r for r in thresholds if r['scenario'] == scenario and r['population'] == population and r['metric'] == metric and float(r['threshold']) == threshold]
                    n = len(tt)
                    assert n == len(cc)
                    st_counts = Counter(r['stochastic_only_status'] for r in tt)
                    prod_counts = Counter(r['stability_class'] for r in tt)
                    survival.append(dict(scenario=scenario, population=population, metric=metric, threshold=threshold, total=n,
                        stochastic_resolved_concordant=st_counts['RESOLVED_CONCORDANT'], stochastic_conflicting_direction=st_counts['CONFLICT'], stochastic_unresolved=st_counts['UNRESOLVED'],
                        stochastic_not_resolved_concordant_proportion=1-st_counts['RESOLVED_CONCORDANT']/n,
                        agreement=prod_counts['AGREEMENT'], resolved_class_disagreement=prod_counts['CHANGE'], unresolved=prod_counts['UNRESOLVED'],
                        unresolved_proportion=prod_counts['UNRESOLVED']/n))
                    assert sum(st_counts.values()) == sum(prod_counts.values()) == n
    write('threshold_survival.csv', survival)
    write('continuous_support_distribution_summary.csv', distributions)
    references = {(r['scenario'], r['geo'], r['year'], r['metric']): r for r in thresholds if float(r['threshold']) == .95}
    support_refs = {(r['scenario'], r['geo'], r['year'], r['metric']): r for r in claims}
    changes = []
    for (scenario, geo, year, metric), r in references.items():
        if scenario != U:
            continue
        b = references[(B, geo, year, metric)]
        su, sb = [support_refs[(sc, geo, year, metric)] for sc in (U, B)]
        changes.append(dict(population=pop(geo), geo=geo, year=year, metric=metric,
            U_source_aware_class=r['source_aware_class'], B_source_aware_class=b['source_aware_class'],
            U_stability_class=r['stability_class'], B_stability_class=b['stability_class'],
            U_stochastic_status=r['stochastic_only_status'], B_stochastic_status=b['stochastic_only_status'],
            source_aware_class_changed=r['source_aware_class'] != b['source_aware_class'],
            final_class_changed=r['stability_class'] != b['stability_class'],
            newly_source_aware_resolved=r['source_aware_class'] == 'UNRESOLVED' and b['source_aware_class'] != 'UNRESOLVED',
            lost_source_aware_resolution=r['source_aware_class'] != 'UNRESOLVED' and b['source_aware_class'] == 'UNRESOLVED',
            newly_final_resolved=r['stability_class'] == 'UNRESOLVED' and b['stability_class'] != 'UNRESOLVED',
            lost_final_resolution=r['stability_class'] != 'UNRESOLVED' and b['stability_class'] == 'UNRESOLVED',
            U_preferred_class_status=su['preferred_class_status'], B_preferred_class_status=sb['preferred_class_status'],
            U_weakest_variant_support=su['weakest_variant_support'], B_weakest_variant_support=sb['weakest_variant_support'],
            paired_support_change=float(sb['weakest_variant_support'])-float(su['weakest_variant_support']) if
                su['weakest_variant_support'] and sb['weakest_variant_support'] else ''))
    write('scenario_classification_changes.csv', changes)
    medium_candidate()
    multispan_summaries()
    verify_serialized()


def multispan(panel, anchors):
    rows = []
    for span in (4, 8):
        half = span // 2
        for geo in P.GEOS:
            for centre in range(2010 + half, 2023 - half, 2):
                years = (centre-half, centre, centre+half)
                if all(P.admissible(anchors, geo, y) for y in years):
                    left, target, right = [P.cmur(panel, geo, y) for y in years]
                    residual = target - (left+right)/2
                    rows.append(dict(population=pop(geo), geo=geo, centre_year=centre, outer_span_years=span,
                        left_year=years[0], right_year=years[2], left_cmur=left, target_cmur=target, right_cmur=right,
                        midpoint_prediction=(left+right)/2, signed_residual=residual, absolute_residual=abs(residual), paired_same_geo_centre=False))
    keys = {span: {(r['geo'], r['centre_year']) for r in rows if r['outer_span_years'] == span} for span in (4, 8)}
    paired = keys[4] & keys[8]
    for r in rows:
        r['paired_same_geo_centre'] = (r['geo'], r['centre_year']) in paired
    write('multispan_pseudogap_diagnostic.csv', rows)
    multispan_summaries()


def multispan_summaries():
    rows = read(HERE / 'multispan_pseudogap_diagnostic.csv')
    summaries = []
    for population in POPS:
        for subset in ('all_eligible', 'paired_same_geo_centre'):
            for span in (4, 8):
                rr = [r for r in rows if r['population'] == population and int(r['outer_span_years']) == span and
                      (subset == 'all_eligible' or r['paired_same_geo_centre'] == 'True')]
                vals = np.array([float(r['signed_residual']) for r in rr])
                absolute = np.abs(vals)
                summaries.append(dict(population=population, subset=subset, outer_span_years=span, eligible_count=len(rr),
                    country_coverage=len({r['geo'] for r in rr}), geographies=';'.join(sorted({r['geo'] for r in rr})),
                    median_absolute_residual=float(np.median(absolute)), mean_absolute_residual=float(np.mean(absolute)),
                    rmse=float(np.sqrt(np.mean(vals**2))), signed_residual_MAD=float(np.median(np.abs(vals-np.median(vals)))),
                    absolute_residual_MAD=float(np.median(np.abs(absolute-np.median(absolute))))))
    write('multispan_summary.csv', summaries)
    paired = defaultdict(dict)
    for r in rows:
        if r['paired_same_geo_centre'] == 'True':
            paired[(r['geo'], r['centre_year'])][r['outer_span_years']] = r
    pairs = []
    for (geo, centre), vv in paired.items():
        absolute = {}
        for span in ('4', '8'):
            r = vv[span]
            # Exact decimal ordering at the released observation precision; no floating tolerance.
            absolute[span] = abs(Decimal(r['target_cmur'])-(Decimal(r['left_cmur'])+Decimal(r['right_cmur']))/2)
        diff = absolute['8']-absolute['4']
        pairs.append(dict(population=pop(geo), geo=geo, centre_year=centre, absolute_residual_4_exact=str(absolute['4']),
            absolute_residual_8_exact=str(absolute['8']), difference_8_minus_4_exact=str(diff),
            ordering='8_GREATER' if diff > 0 else '4_GREATER' if diff < 0 else 'EQUAL'))
    write('multispan_paired_comparison.csv', pairs)


def medium_candidate():
    prior = read(PRIOR / 'medium_term_harmonized_rules.csv')
    rows = []
    for rule, role in (('COMMON_TWO_SIDED_95', 'candidate_PRIMARY'), ('COMMON_DIRECTIONAL_95', 'SECONDARY_robustness')):
        for population in POPS:
            rr = [r for r in prior if r['rule'] == rule and r['population'] == population]
            wa = np.array([float(r['A_ci_high'])-float(r['A_ci_low']) for r in rr])
            wb = np.array([float(r['B1_ci_high'])-float(r['B1_ci_low']) for r in rr])
            rows.append(dict(role=role, rule=rule, population=population,
                agreement=sum(r['stability_class'] == 'AGREEMENT' for r in rr), resolved_class_disagreement=sum(r['stability_class'] == 'CHANGE' for r in rr),
                unresolved=sum(r['stability_class'] == 'UNRESOLVED' for r in rr),
                A_unresolved=sum(r['A_class'] == 'UNRESOLVED' for r in rr), B1_unresolved=sum(r['B1_class'] == 'UNRESOLVED' for r in rr),
                median_A_CI_width=float(np.median(wa)), median_B1_CI_width=float(np.median(wb)), median_paired_B1_A_CI_width_ratio=float(np.median(wb/wa)),
                A_changed_geographies=';'.join(r['geo'] for r in rr if r['A_class_changed'] == 'True'),
                B1_changed_geographies=';'.join(r['geo'] for r in rr if r['B1_class_changed'] == 'True'),
                comparison_changed_geographies=';'.join(r['geo'] for r in rr if r['stability_class_changed'] == 'True')))
    write('medium_term_revision_candidate.csv', rows)


def verify_serialized():
    continuous = read(HERE / 'calibration_scenario_support.csv')
    cc = read(HERE / 'claim_level_continuous_support.csv')
    thresholds = read(HERE / 'threshold_claim_classifications.csv')
    survival = read(HERE / 'threshold_survival.csv')
    changes = read(HERE / 'scenario_classification_changes.csv')
    lookup = {(r['scenario'], r['geo'], r['year'], r['metric'], r['variant']): r for r in continuous}
    assert len(lookup) == len(continuous) == 7616
    assert len(cc) == len({(r['scenario'], r['geo'], r['year'], r['metric']) for r in cc}) == 1904
    assert len(thresholds) == len({(r['scenario'], r['geo'], r['year'], r['metric'], r['threshold']) for r in thresholds}) == 13328
    assert len(changes) == len({(r['geo'], r['year'], r['metric']) for r in changes}) == 952
    for r in continuous:
        if r['metric'] == 'turning' and r['stochastic'] == 'True' and r['estimable'] == 'True':
            for prefix, year in (('earlier', int(r['year'])-1), ('later', int(r['year']))):
                annual = lookup[(r['scenario'], r['geo'], str(year), 'annual_change', r['variant'])]
                for field in ('p_positive', 'p_negative'):
                    assert float(r[prefix+'_change_'+field]) == float(annual[field])
                assert float(r[prefix+'_change_support']) == float(annual['preferred_sign_support'])
                assert r[prefix+'_change_preferred_sign'] == annual['preferred_sign']
                assert r[prefix+'_change_class_0_95'] == annual['variant_class']
    frozen = {(r['geo'], r['year'], r['metric']): r for r in read(P.OUT / 'temporal_stability_comparison.csv') if r['metric'] in METRICS and int(r['year']) <= 2022}
    for r in thresholds:
        votes = []
        stoch = []
        for v in P.VARIANTS:
            source = lookup[(r['scenario'], r['geo'], r['year'], r['metric'], v)]
            # Independent threshold implementation uses preferred labels plus continuous support;
            # turning minimum of constituent supports is sufficient only alongside its preferred class.
            if source['variant_class'] == 'NOT_APPLICABLE':
                label = 'NOT_APPLICABLE'
            elif source['source_aware_estimate'] == '':
                label = 'UNRESOLVED'
            elif v == 'deterministic_midpoint':
                label = source['variant_class']
            else:
                label = source['preferred_class'] if source['preferred_class'] != 'TIE' and float(source['preferred_sign_support']) >= float(r['threshold']) else 'UNRESOLVED'
            assert label == r[v+'_class']
            if label != 'NOT_APPLICABLE':
                votes.append(label)
                if v in STOCH:
                    stoch.append(label)
        resolved_labels = {x for x in stoch if x != 'UNRESOLVED'}
        status = 'CONFLICT' if len(resolved_labels) > 1 else 'RESOLVED_CONCORDANT' if len(stoch) >= 2 and all(x != 'UNRESOLVED' for x in stoch) else 'UNRESOLVED'
        final = votes[0] if len(votes) >= 3 and len(set(votes)) == 1 and votes[0] != 'UNRESOLVED' else 'UNRESOLVED'
        assert status == r['stochastic_only_status'] and final == r['source_aware_class']
        assert compare(r['branchA_class'], final) == r['stability_class']
        if r['scenario'] == U and float(r['threshold']) == .95:
            f = frozen[(r['geo'], r['year'], r['metric'])]
            assert (r['source_aware_class'], r['stability_class']) == (f['source_aware_class'], f['stability_class'])
    for r in cc:
        ss = [lookup[(r['scenario'], r['geo'], r['year'], r['metric'], v)] for v in STOCH]
        ss = [s for s in ss if s['applicable'] == 'True']
        labels = {s['preferred_class'] for s in ss}
        status = 'MISSING_COVERAGE' if any(s['estimable'] != 'True' for s in ss) else 'PREFERRED_TIE' if 'TIE' in labels else 'CONCORDANT' if len(labels) == 1 else 'CONFLICT'
        assert status == r['preferred_class_status']
        if status == 'CONCORDANT':
            assert float(r['weakest_variant_support']) == min(float(s['preferred_sign_support']) for s in ss)
        else:
            assert not r['weakest_variant_support']
    for row in survival:
        rr = [r for r in thresholds if all(r[k] == row[k] for k in ('scenario', 'population', 'metric', 'threshold'))]
        assert int(row['total']) == len(rr)
        for name, status in (('stochastic_resolved_concordant', 'RESOLVED_CONCORDANT'), ('stochastic_conflicting_direction', 'CONFLICT'), ('stochastic_unresolved', 'UNRESOLVED')):
            assert int(row[name]) == sum(r['stochastic_only_status'] == status for r in rr)
        for name, status in (('agreement', 'AGREEMENT'), ('resolved_class_disagreement', 'CHANGE'), ('unresolved', 'UNRESOLVED')):
            assert int(row[name]) == sum(r['stability_class'] == status for r in rr)
        assert int(row['stochastic_resolved_concordant']) == int(row['agreement']) + int(row['resolved_class_disagreement'])
    common_keys = {(r['geo'], r['year'], r['metric']) for r in cc if r['scenario'] == U and r['preferred_class_status'] == 'CONCORDANT'} & {
                   (r['geo'], r['year'], r['metric']) for r in cc if r['scenario'] == B and r['preferred_class_status'] == 'CONCORDANT'}
    for row in read(HERE / 'continuous_support_distribution_summary.csv'):
        selected = [r for r in cc if all(r[k] == row[k] for k in ('scenario', 'population', 'metric')) and r['weakest_variant_support'] != '' and
                    (row['subset'] == 'scenario_concordant' or (r['geo'], r['year'], r['metric']) in common_keys)]
        vals = np.array([float(r['weakest_variant_support']) for r in selected])
        assert len(vals) == int(row['concordant_with_support'])
        assert list(np.quantile(vals, [.1, .25, .5, .75, .9])) == [float(row[k]) for k in ('p10', 'q25', 'median', 'q75', 'p90')]
        for cut in (.8, .9, .95, .99):
            assert float(row[f'proportion_ge_{cut:.2f}']) == np.mean(vals >= cut)
    for population in POPS:
        for metric in METRICS:
            rr = [r for r in changes if r['population'] == population and r['metric'] == metric]
            u = next(r for r in survival if r['scenario'] == U and r['population'] == population and r['metric'] == metric and float(r['threshold']) == .95)
            b = next(r for r in survival if r['scenario'] == B and r['population'] == population and r['metric'] == metric and float(r['threshold']) == .95)
            new = sum(r['newly_final_resolved'] == 'True' for r in rr)
            lost = sum(r['lost_final_resolution'] == 'True' for r in rr)
            assert new-lost == int(u['unresolved'])-int(b['unresolved'])
            assert sum(r['U_stability_class'] == 'UNRESOLVED' for r in rr) == int(u['unresolved'])
            assert sum(r['B_stability_class'] == 'UNRESOLVED' for r in rr) == int(b['unresolved'])
    audit = read(HERE / 'paired_draw_audit.csv')
    assert all(r['U_identity_sha256'] == r['B_identity_sha256'] and r['U_sign_sha256'] == r['B_sign_sha256'] and
               r['exact_amplitude_transform'] == 'True' and r['independent_stream_replay_equal'] == 'True' for r in audit)
    assert sha(HERE / 'scenario_u_production_reproduction.csv') == sha(P.OUT / 'b2r_short_window_support.csv')
    verify_paired_draws(lookup)
    verify_multispan()
    verify_prior_details()
    before = json.loads((HERE / 'input_integrity.json').read_text())['before_sha256']
    after = inventory()
    changed = [p for p in before if before[p] != after.get(p)]
    added = sorted(set(after)-set(before))
    assert not changed and not added
    save('input_integrity.json', dict(before_sha256=before, after_sha256=after, changed=changed, added_outside_diagnostic_directory=added,
        all_preexisting_files_unchanged=True, preexisting_files_checked=len(before), git_available=False))
    v = json.loads((HERE / 'verification.json').read_text())
    v.update(serialized_row_counts=dict(variant_support=len(continuous), claim_support=len(cc), threshold_claims=len(thresholds), scenario_pairs=len(changes)),
        all_threshold_classes_independently_verified=True, all_survival_counts_reconcile=True, scenario_changes_reconcile=True,
        paired_draw_hashes_independently_replayed=True, B_sign_support_independently_recomputed=True,
        multispan_admissibility_and_residuals_verified=True, continuous_distribution_statistics_verified=True,
        all_prior_P1A_headline_aggregates_reproduced=True, turning_constituents_match_annual_support_rows=True,
        conditional_support_and_conflicts_verified=True, original_files_checked=len(before), original_files_changed=0,
        scenario_U_production_sha256=sha(P.OUT / 'b2r_short_window_support.csv'), pass_all=True)
    save('verification.json', v)
    print('PASS serialized row verification, all threshold aggregates, paired scenarios and original-file integrity.', flush=True)


def verify_paired_draws(lookup):
    """Replay seeds afresh and check saved draw hashes, B supports and unchanged anchors."""
    panel, anchors = P.load_panel()
    _, by_geo, pooled, eligible = P.prepare(panel, anchors)
    audit = {(r['geo'], r['variant'], int(r['year'])): r for r in read(HERE / 'paired_draw_audit.csv')}
    for geo in P.GEOS:
        for variant in P.VARIANTS:
            original = P.simulate(panel, anchors, by_geo, pooled, eligible, geo, variant)
            if original is None:
                continue
            det = variant == 'deterministic_midpoint'
            n = 1 if det else P.N_DRAWS
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([P.SEED, P.GEOS.index(geo), P.VARIANTS.index(variant)])))
            scaled_path = {}
            for year in P.YEARS:
                if P.admissible(anchors, geo, year):
                    scaled_path[year] = np.full(n, P.cmur(panel, geo, year), dtype=float)
                    assert np.array_equal(scaled_path[year], original[year])
                elif (geo, year) in eligible:
                    if det:
                        ids, signs, raw = np.zeros(n, dtype=np.int64), np.ones(n), np.zeros(n)
                    else:
                        pool = np.asarray(by_geo[geo], dtype=float) if variant == 'country_signed' else pooled
                        if variant == 'pooled_absolute_random_sign':
                            pool = np.abs(pool)
                        ids = rng.choice(len(pool), size=n, replace=True)
                        signs = rng.choice(np.array([-1., 1.]), size=n, replace=True) if variant == 'pooled_absolute_random_sign' else np.ones(n)
                        raw = pool[ids] * signs
                    record = audit[(geo, variant, year)]
                    assert array_hash(ids) == record['U_identity_sha256'] == record['B_identity_sha256']
                    assert array_hash(signs) == record['U_sign_sha256'] == record['B_sign_sha256']
                    scaled = raw * float(1 / np.sqrt(2))
                    assert array_hash(raw) == record['U_residual_sha256']
                    assert array_hash(scaled) == record['B_residual_sha256']
                    midpoint = eligible[(geo, year)]
                    assert np.array_equal(np.clip(midpoint + raw, 0, 100), original[year])
                    scaled_path[year] = np.clip(midpoint + scaled, 0, 100)
                    assert int(record['U_clipped_draws']) == np.count_nonzero((midpoint+raw < 0) | (midpoint+raw > 100))
                    assert int(record['B_clipped_draws']) == np.count_nonzero((midpoint+scaled < 0) | (midpoint+scaled > 100))
            for metric in METRICS:
                start = 2011 if metric == 'annual_change' else 2012
                for year in range(start, 2023):
                    saved = lookup[(B, geo, str(year), metric, variant)]
                    required = (year-1, year) if metric == 'annual_change' else (year-2, year-1, year)
                    if not all(y in scaled_path for y in required):
                        assert saved['source_aware_estimate'] == ''
                        continue
                    if metric == 'annual_change':
                        values = scaled_path[year]-scaled_path[year-1]
                    elif metric == 'acceleration':
                        values = np.zeros(1) if det and (geo, year-1) in eligible else scaled_path[year]-2*scaled_path[year-1]+scaled_path[year-2]
                    else:
                        values = ((scaled_path[year]-scaled_path[year-1])*(scaled_path[year-1]-scaled_path[year-2]) < 0).astype(float)
                    assert float(saved['source_aware_estimate']) == float(np.mean(values))
                    if not det and metric != 'turning':
                        assert float(saved['p_positive']) == np.count_nonzero(values > 0)/n
                        assert float(saved['p_negative']) == np.count_nonzero(values < 0)/n
                    if det:
                        unscaled = lookup[(U, geo, str(year), metric, variant)]
                        assert all(saved[k] == unscaled[k] for k in saved if k not in ('scenario', 'residual_amplitude_scale'))


def verify_multispan():
    panel, anchors = P.load_panel()
    rows = read(HERE / 'multispan_pseudogap_diagnostic.csv')
    expected = set()
    for span in (4, 8):
        half = span//2
        for geo in P.GEOS:
            for year in range(2010+half, 2023-half, 2):
                if all(P.admissible(anchors, geo, t) for t in (year-half, year, year+half)):
                    expected.add((geo, year, span))
    assert len(rows) == len(expected) == len({(r['geo'], int(r['centre_year']), int(r['outer_span_years'])) for r in rows})
    assert expected == {(r['geo'], int(r['centre_year']), int(r['outer_span_years'])) for r in rows}
    original = {(r['geo'], r['interior_anchor_year']): r for r in read(P.OUT / 'b2r_residual_pool.csv')}
    for r in rows:
        years = [int(r[k]) for k in ('left_year', 'centre_year', 'right_year')]
        assert all(P.admissible(anchors, r['geo'], year) for year in years)
        v = [P.cmur(panel, r['geo'], year) for year in years]
        residual = v[1]-(v[0]+v[2])/2
        assert residual == float(r['signed_residual'])
        if r['outer_span_years'] == '4':
            assert residual == float(original[(r['geo'], r['centre_year'])]['signed_residual_pp'])
    for r in read(HERE / 'multispan_summary.csv'):
        selected = [x for x in rows if x['population'] == r['population'] and x['outer_span_years'] == r['outer_span_years'] and
                    (r['subset'] == 'all_eligible' or x['paired_same_geo_centre'] == 'True')]
        vals = np.array([float(x['signed_residual']) for x in selected])
        assert len(vals) == int(r['eligible_count'])
        assert float(r['median_absolute_residual']) == np.median(np.abs(vals))
        assert float(r['mean_absolute_residual']) == np.mean(np.abs(vals))
        assert float(r['rmse']) == np.sqrt(np.mean(vals**2))


def verify_prior_details():
    """Read-only replay of all P1A influence rows and aggregate cause/influence headlines."""
    summary = json.loads((PRIOR / 'aggregate_summary.json').read_text())
    causes = read(PRIOR / 'short_window_unresolved_causes.csv')
    influence = read(PRIOR / 'short_window_variant_influence.csv')
    frozen = defaultdict(dict)
    for r in read(P.OUT / 'b2r_short_window_support.csv'):
        frozen[(r['geo'], r['year'], r['metric'])][r['variant']] = r
    for r in influence:
        vv = frozen[(r['geo'], r['year'], r['metric'])]
        labels = [x['variant_class'] for v, x in vv.items() if v != r['omitted_voting_variant'] and x['variant_class'] != 'NOT_APPLICABLE']
        cls = labels[0] if len(labels) >= int(r['minimum_applicable_variants']) and len(set(labels)) == 1 and labels[0] != 'UNRESOLVED' else 'UNRESOLVED'
        assert cls == r['source_aware_class']
        assert compare(r['branchA_class'], cls) == r['stability_class']
    for row in summary['influence']:
        rr = [r for r in influence if all(r[k] == row[k] for k in ('population', 'period', 'metric', 'diagnostic'))]
        for label in ('AGREEMENT', 'CHANGE', 'UNRESOLVED'):
            assert row[label] == sum(r['stability_class'] == label for r in rr)
        for output, field in (('source_aware_changed', 'source_aware_class_changed'), ('comparison_changed', 'stability_class_changed'),
                              ('newly_source_aware_resolved', 'newly_source_aware_resolved'), ('newly_comparison_resolved', 'newly_comparison_resolved')):
            assert row[output] == sum(r[field] == 'True' for r in rr)
    for r in causes:
        vv = frozen[(r['geo'], r['year'], r['metric'])]
        for v in P.VARIANTS:
            assert r[v+'_class'] == vv[v]['variant_class'] and r[v+'_reason'] == vv[v]['reason']
    for row in summary['causes']:
        rr = [r for r in causes if all(r[k] == row[k] for k in ('population', 'period', 'metric'))]
        assert row['unresolved'] == len(rr)
        assert row['primary'] == dict(Counter(r['primary_cause'] for r in rr))
        assert row['flags'] == {k: sum(r[k] == 'True' for r in rr) for k in row['flags']}


def md_table(headers, rows):
    escape = lambda value: str(value).replace('|', '\\|')
    return '\n'.join(['| ' + ' | '.join(map(escape, headers)) + ' |', '| ' + ' | '.join(['---']*len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(escape, r)) + ' |' for r in rows])


def report():
    survival = read(HERE / 'threshold_survival.csv')
    distributions = read(HERE / 'continuous_support_distribution_summary.csv')
    changes = read(HERE / 'scenario_classification_changes.csv')
    claims = read(HERE / 'claim_level_continuous_support.csv')
    thresholds = read(HERE / 'threshold_claim_classifications.csv')
    spans = read(HERE / 'multispan_summary.csv')
    pairs = read(HERE / 'multispan_paired_comparison.csv')
    medium = read(HERE / 'medium_term_revision_candidate.csv')
    verification = json.loads((HERE / 'verification.json').read_text())
    assert verification['pass_all']
    short_names = {U: 'U: unscaled stress', B: 'B: Brownian SD'}
    def selected(scenario, population, metric, threshold):
        return next(r for r in survival if r['scenario'] == scenario and r['population'] == population and
                    r['metric'] == metric and float(r['threshold']) == threshold)
    def adu(r):
        return '/'.join(r[k] for k in ('agreement', 'resolved_class_disagreement', 'unresolved'))
    reference_rows = []
    for population in POPS:
        for metric in METRICS:
            u, b = [selected(s, population, metric, .95) for s in SCALES]
            delta = 100*(float(b['unresolved_proportion'])-float(u['unresolved_proportion']))
            changed = [r for r in changes if r['population'] == population and r['metric'] == metric]
            newly = sum(r['newly_final_resolved'] == 'True' for r in changed)
            lost = sum(r['lost_final_resolution'] == 'True' for r in changed)
            reference_rows.append([population, metric, adu(u), adu(b), newly, lost, f'{delta:+.2f}'])
    reference_table = md_table(['Population', 'Metric', 'U A/D/U', 'B A/D/U', 'New resolutions', 'Lost', 'Change unresolved, pp'], reference_rows)
    range_rows = []
    for scenario in SCALES:
        for metric in METRICS:
            rr = [r for r in survival if r['scenario'] == scenario and r['population'] == POPS[0] and r['metric'] == metric]
            assert all(float(r['stochastic_not_resolved_concordant_proportion']) > .5 for r in rr)
            assert all(float(r['unresolved_proportion']) > .5 for r in rr)
            range_rows.append([short_names[scenario], metric,
                f"{100*min(float(r['unresolved_proportion']) for r in rr):.1f}–{100*max(float(r['unresolved_proportion']) for r in rr):.1f}%",
                selected(scenario, POPS[0], metric, .80)['stochastic_resolved_concordant'],
                selected(scenario, POPS[0], metric, .99)['stochastic_resolved_concordant']])
    range_table = md_table(['Scenario', 'Metric', 'Final unresolved range over grid', 'Stochastic concordant at .80', 'At .99'], range_rows)
    dist_table = md_table(['Scenario', 'Metric', 'n concordant', 'Preferred conflicts', 'P10', 'Q25', 'Median', 'Q75', 'P90', '≥.80', '≥.90', '≥.95', '≥.99'],
        [[short_names[r['scenario']], r['metric'], r['concordant_with_support'], r['conflict']] +
         [f"{float(r[k]):.3f}" for k in ('p10', 'q25', 'median', 'q75', 'p90')] +
         [f"{100*float(r[f'proportion_ge_{cut:.2f}']):.1f}%" for cut in (.8, .9, .95, .99)]
         for r in distributions if r['population'] == POPS[0] and r['subset'] == 'scenario_concordant'])
    span_table = md_table(['Population', 'Subset', 'Span', 'n', 'Geographies', 'Median |e|', 'Mean |e|', 'RMSE', 'Signed MAD'],
        [[r['population'], r['subset'], r['outer_span_years'], r['eligible_count'], r['country_coverage']] +
         [f"{float(r[k]):.3f}" for k in ('median_absolute_residual', 'mean_absolute_residual', 'rmse', 'signed_residual_MAD')] for r in spans])
    medium_table = md_table(['Candidate role', 'Population', 'A/D/U', 'A/B1 unresolved', 'Median A/B1 CI width', 'Median paired ratio', 'Changed geographies'],
        [[r['role'], r['population'], adu(r), r['A_unresolved']+'/'+r['B1_unresolved'],
          f"{float(r['median_A_CI_width']):.4f}/{float(r['median_B1_CI_width']):.4f}",
          f"{float(r['median_paired_B1_A_CI_width_ratio']):.3f}",
          'A: '+(r['A_changed_geographies'] or 'none')+'; B1: '+(r['B1_changed_geographies'] or 'none')] for r in medium])
    ms_pairs = [r for r in pairs if r['population'] == POPS[0]]
    order = Counter(r['ordering'] for r in ms_pairs)
    disagreement_rows = [r for r in thresholds if r['stability_class'] == 'CHANGE']
    disagreements = '; '.join(f"{geo} {year} {metric}, thresholds " + ', '.join(sorted({r['threshold'] for r in disagreement_rows if (r['geo'],r['year'],r['metric']) == (geo,year,metric)}, key=float))
        for geo, year, metric in sorted({(r['geo'],r['year'],r['metric']) for r in disagreement_rows})) or 'none'
    changed_count = sum(r['final_class_changed'] == 'True' for r in changes)
    lost_count = sum(r['lost_final_resolution'] == 'True' for r in changes)
    common_text = '; '.join(f"{metric}: " + ' → '.join(f"{float(next(r for r in distributions if r['scenario'] == scenario and r['population'] == POPS[0] and r['metric'] == metric and r['subset'] == 'common_concordant_pairs')['median']):.3f}" for scenario in SCALES) for metric in METRICS)
    eu_medians = '; '.join(f"{metric}: " + ' → '.join(f"{float(next(r for r in distributions if r['scenario'] == scenario and r['population'] == POPS[1] and r['metric'] == metric and r['subset'] == 'scenario_concordant')['median']):.3f}" for scenario in SCALES) for metric in METRICS)
    text = f'''# Scientific decision evidence: calibration robustness

10 November 2026. U is the PRIMARY empirical unscaled four-year stress test. B is the SECONDARY Brownian-model amplitude scenario, multiplying the same sampled residuals by 1/√2. Neither scenario supplies calibrated two-year uncertainty probabilities. A/D/U means agreement / resolved-class disagreement / unresolved against the released series. This report records diagnostics; it makes no manuscript, publication, novelty, P4 or journal decision.

## 1. Does the short-window boundary survive Brownian scaling?

Yes, within these two scenarios and the pre-specified threshold grid. Most claims lack resolved cross-variant sign support even after amplitude reduction. The statement concerns CMUR under the tested calibration families and consensus design, not the unknown distribution of true annual errors.

## 2. What changes at the 0.95 reference threshold?

{reference_table}

All {changed_count} changed final comparisons become agreements; {lost_count} lose resolution. Exact geography/year/metric identities and old/new classes are in scenario_classification_changes.csv (filter final_class_changed=True). Scaling changes a meaningful set of individual classifications but reduces unresolvedness by only the percentage-point amounts above. Annual change is most sensitive in absolute counts and percentage-point change; acceleration is least sensitive by those measures. Relative growth from tiny baseline agreement counts gives a different ranking and should not be confused with panel-level sensitivity. EU27_2020 gains no resolution.

## 3. Is the conclusion dominated by choosing 0.95?

No majority conclusion depends uniquely on 0.95: every metric remains more than half unresolved throughout 0.80, 0.85, 0.90, 0.925, 0.95, 0.975 and 0.99, in both scenarios. Detailed counts change with threshold:

{range_table}

All threshold-level stochastic conflict counts are {max(int(r['stochastic_conflicting_direction']) for r in survival)}: the opposing preferred classes described below do not both reach these thresholds. Stochastic-only concordance and production-style final resolution counts coincide in these simulations across this grid; keeping the midpoint vote adds no veto at these thresholds. The grid is a robustness diagnosis, not a threshold-selection exercise. Resolved disagreement does occur under B at low thresholds: {disagreements}. This annual change is conventional negative versus source-aware positive; there is no final resolved disagreement at 0.95. EU27_2020 remains entirely unresolved throughout the grid in both scenarios.

## 4–5. What do continuous supports reveal?

The weakest-variant statistic is the minimum preferred-sign support across applicable stochastic variants only when their preferred classes agree. It is a robustness statistic, not a joint or calibrated probability. For turning it is the minimum of constituent-change supports across variants, with the preferred TURN/NO_TURN class inferred from the constituent preferred signs. Strict reversal-event frequencies are stored separately and never substituted for this support.

Member-State distributions conditional on concordant preferred classes are:

{dist_table}

The distribution denominators differ because scaling can change preferred-class concordance. Missing-coverage counts are unchanged ({'/'.join(str(next(r for r in distributions if r['scenario'] == U and r['population'] == POPS[0] and r['metric'] == metric and r['subset'] == 'scenario_concordant')['missing_coverage']) for metric in METRICS)} by metric); no preferred ties occur. Conflicting preferred classes are explicitly excluded from a single weakest-support number, not averaged. The common-concordant paired subset is also saved: its medians U → B are {common_text}. EU27_2020 conditional medians U → B are {eu_medians}; none reaches cross-variant support of 0.80. These graded shifts and numerous weak direction conflicts are obscured by a single 0.95 class. Acceleration shows especially low support and frequent preferred-direction conflicts.

## 6. What does the exploratory eight-year exercise show?

Task F is feasible for a descriptive Member-State exercise. All predictions use genuinely admissible outer and target anchors. Values below are in CMUR percentage points; MAD is the unscaled median absolute deviation about the median signed residual.

{span_table}

The paired Member-State subset has the same geography and central target for both spans, with all five required anchors admissible. Eight-year absolute residuals are larger in {order['8_GREATER']} cases, smaller in {order['4_GREATER']}, and equal in {order['EQUAL']} (ordering evaluated at exact released decimal precision). Thus pooled and paired aggregate scales are larger for eight years, consistent with smaller typical error at shorter observed spans, but casewise monotonicity fails. Nested windows, overlapping targets and within-country dependence prevent treating these as independent replications. EU27_2020 has only {sum(r['population'] == POPS[1] for r in pairs)} pairs and is descriptive only. This comparison does not observe two-year validation errors, validate Brownian dynamics or the 1/√2 factor, calibrate a two-year scale, or support a validated scaling exponent.

## 7. What level of scientific interpretation is supported?

The evidence goes beyond a threshold-specific result: the majority boundary survives the full threshold grid and the approved model-motivated amplitude reduction, with continuous summaries showing predominantly weak support. It supports a **robust CMUR-specific temporal-support boundary within this tested sensitivity design**. The computational exercise also illustrates a method for separating coverage, preferred-direction conflicts, support thresholds and amplitude assumptions; it does not establish a general empirical boundary for other indicators. Calibration from four-year residuals to actual two-year error remains unidentified without further assumptions or validation. The Controller retains all final contribution, novelty and P4 decisions.

## Medium-term revision candidate (no new models)

{medium_table}

COMMON_TWO_SIDED_95 is the candidate symmetric primary; COMMON_DIRECTIONAL_95 is secondary robustness. BE changes B1 positive → unresolved; FI and LU negative → unresolved under the candidate primary, changing their comparisons to unresolved. SI changes A unresolved → positive under the secondary, while its comparison remains unresolved. These rows are reused from verified P1A fits; no medium-term model was refitted here.

## Verification

Production and P1A gates passed before alternatives. Scenario U's complete production support CSV is byte-identical (including endpoints and all stochastic fields/classes); B alone changes residual amplitude. Independent seed replay verifies residual identities/signs/order and every stored B sign support. All {len(thresholds)} threshold claim rows are independently reconstructed from serialized continuous support, and all reported aggregate tables are generated from serialized files. Existing verification suites passed in isolation; {verification['original_files_checked']} pre-existing files, including all P1A outputs, manuscript, figures, frozen sources and Zenodo evidence, remain byte-identical. No scientific or technical blocker remains; the workspace lacks Git metadata, so integrity is verified by SHA-256 rather than a Git diff. Existing integration checks use the preserved build logs/PDF, not a new LaTeX build.

The most decision-relevant finding is that Brownian scaling raises 0.95 resolutions but leaves the majority short-window boundary intact even at 0.80, so the qualitative finding is not an artifact of the 0.95 cutoff alone.
'''
    (HERE / 'SCIENTIFIC_DECISION_SUMMARY.md').write_text(text, encoding='utf-8')
    figures()
    files = sorted(p.name for p in HERE.iterdir() if p.is_file())
    if 'README.md' not in files:
        files.append('README.md')
        files.sort()
    readme = f'''# WP P1C: CMUR calibration robustness

Completed 10 November 2026. This directory contains diagnostics only. No existing manuscript, data, CSV, figure, P1A result or Zenodo evidence was edited. No external data or literature was used. The uncertainty interpretation follows the Controller's fixed scientific basis.

## Authority and gates

Production authority is scripts/wp1b_r.py, whose OUT, load_panel, prepare, simulate, classify_array, short_metric and conventional_class functions locate and define the frozen inputs and mechanics. Its data/processed inputs are cmur_source_provenance_panel.csv, verified_anchor_map.csv, b2r_residual_pool.csv, b2r_short_window_support.csv and temporal_stability_comparison.csv. Medium-term rows come from the prior harmonized CSV, checked against the frozen branchA_medium_term.csv and branchB1_anchor_medium_term.csv located by production. The manuscript pipeline and scripts/verify_wp1a.py, scripts/verify_wp1b_r.py and manuscript/scripts/verify_integration.py establish production reproduction. No medium-term model was refitted.

Every mandatory P1A input was read: baseline_reproduction.csv, medium_term_harmonized_rules.csv, short_window_unresolved_causes.csv, short_window_variant_influence.csv, SCIENTIFIC_DIAGNOSTIC_SUMMARY.md, aggregate_summary.json, verification.json and independent_verification.json. Production gates reproduce 15/0/12 medium-term, 10/0/314 annual change, 2/0/295 turning, and 1/0/296 acceleration for Member States. P1A symmetric rules reproduce 12/0/15 and 15/0/12; both strict and exploratory midpoint removals reproduce zero new acceleration resolutions, for Member States and EU separately.

## Scenarios and paired simulation

- UNSCALED_STRESS (U) is PRIMARY: the empirical four-year pseudo-gap residual stress-test amplitude, with no claim to calibrated two-year errors or formal upper confidence bounds.
- BROWNIAN_SD_1_SQRT2 (B) is SECONDARY: exactly float(1/np.sqrt(2)) = {SCALES[B]!r} multiplies each sampled residual before adding its midpoint and clipping. This is the Controller-specified Brownian-motion/bridge SD scenario, not the true CMUR correction.

The fixed seed is {P.SEED}, with {P.N_DRAWS} trajectories and the unchanged PCG64 SeedSequence([seed, geography index, variant index]) streams. Residual pool order, sampled identities, randomized signs, draw order, independent odd-year draws, clipping to [0,100], anchor values and country applicability are unchanged. A single identity/sign sequence feeds both scenarios; this also preserves draws with zero-valued or repeated residuals. Integer-index sampling was verified against production's value-sampling API, and every U path was array-equal to production simulate. The full regenerated U CSV is byte-identical to frozen production, not only its classes. Scenario B is evaluated only after this gate passes. No residual ×1/2 scenario was run.

paired_draw_audit.csv stores one row per geography/variant/eligible odd-year block, including original stream indices, hashes of identities, signs and unscaled/scaled pre-clipping residuals, and clipped-draw counts. Identical identity/sign hashes are deliberate: both scenarios use the same sampled objects. Verification independently replays the seeds, reconstructs these hashes and both paths, and checks B estimates and sign frequencies. No random stream is reseeded or rerouted after variant exclusion. Full trajectories are reproducible from the saved code and frozen inputs; they are not duplicated as a large draw archive.

The deterministic midpoint has one path, unchanged in B, and uses the production exact-zero acceleration representation for anchor–midpoint–anchor triples. It is structural information without invented probability values.

## Continuous support definitions

calibration_scenario_support.csv has 7,616 rows: two scenarios × 28 geographies × (12 annual + 11 turning + 11 acceleration) × four variants, through 2022 only. Applicable missing-coverage rows remain in the grid as unresolved; country_signed is NOT_APPLICABLE for CZ and EL. No endpoint claim is silently dropped from the production byte-identity check, but 2023–2024 are outside this WP's scenario/threshold analysis.

For stochastic annual change and acceleration, p_positive and p_negative are the empirical fractions strictly above/below zero; zeros support neither. preferred_sign/class selects the larger fraction, with TIE explicit. preferred_sign_support is max(p_positive,p_negative). distance_from_0_95_threshold is the signed difference support−0.95. Fractions describe the sensitivity ensemble; they are neither posterior probabilities nor calibrated frequentist confidence probabilities.

Turning is not a sign-valued metric: p_positive/p_negative and preferred_sign are intentionally blank. preferred_class is TURN/NO_TURN from the two constituent preferred signs (or TIE). Its preferred_sign_support is the MINIMUM of the two constituent max-sign supports, as indicated by support_definition. The earlier/later p_positive, p_negative, preferred sign, support and 0.95 class are exposed. At each threshold both constituent change classes are recomputed; unresolved either side means turning unresolved. source_aware_estimate retains the production strict reversal-event frequency separately. It is never relabelled as a turning sign-support probability. Deterministic rows have all probability/support quantities blank.

claim_level_continuous_support.csv has 1,904 scenario/claim rows. All applicable stochastic variants are considered, including applicable unresolved ones; NOT_APPLICABLE is omitted. Preferred-class status is MISSING_COVERAGE, PREFERRED_TIE, CONCORDANT or CONFLICT. Only CONCORDANT rows get weakest_variant_support = minimum preferred support across applicable variants. Conflicts receive no averaged probability. Turning concordance concerns TURN/NO_TURN, not equality of both underlying sign directions. This statistic does not assume variants are independent. The deterministic estimate/class is included solely as structural reference.

## Threshold and comparison conventions

The exact grid is {', '.join(map(str, GRID))}. threshold_claim_classifications.csv contains 13,328 scenario/claim/threshold rows, with all four recomputed variant classes. Stochastic-only resolution requires at least TWO applicable stochastic variants, every one resolved and all the same class; normally three vote, but CZ and EL have two. A conflicting-direction claim has at least two incompatible resolved stochastic classes, even if another applicable vote is unresolved; conflict takes precedence in the exclusive status. All other claims are stochastic UNRESOLVED. Preferred conflicts below threshold remain distinct in the continuous output.

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

input_integrity.json records all {verification['original_files_checked']} pre-existing file hashes before and after, including the entire previous diagnostic directory; zero changes and zero additions outside this new directory are required. Verification never writes to a pre-existing path. The supplied workspace has no .git metadata, so Git diff/commit is unavailable; SHA-256 provides whole-file integrity evidence. All deliverables are saved locally. No blocker remains.

## Exact file inventory

''' + '\n'.join('- ' + name for name in files) + '\n'
    (HERE / 'README.md').write_text(readme, encoding='utf-8')
    verify_serialized()
    print('Regenerated scientific report, README and diagnostic figures from serialized outputs.', flush=True)


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False,
        'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.dpi': 300})
    colors = {U: '#1764A2', B: '#C56B16'}
    names = {U: 'U: unscaled stress', B: 'B: Brownian SD'}
    labels = {'annual_change': 'Annual change', 'turning': 'Turning', 'acceleration': 'Acceleration'}
    survival = read(HERE / 'threshold_survival.csv')
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.5), sharey=True, layout='constrained')
    for ax, metric in zip(axes, METRICS):
        for scenario in SCALES:
            rr = sorted([r for r in survival if r['population'] == POPS[0] and r['scenario'] == scenario and r['metric'] == metric], key=lambda r: float(r['threshold']))
            x = [float(r['threshold']) for r in rr]
            ax.plot(x, [100*float(r['stochastic_not_resolved_concordant_proportion']) for r in rr], color=colors[scenario], lw=2, marker='o', ms=3, label=names[scenario])
            ax.plot(x, [100*float(r['unresolved_proportion']) for r in rr], color=colors[scenario], lw=1.6, ls=':', label='_nolegend_')
        ax.set_title(labels[metric], fontweight='bold')
        ax.set_xlabel('Sign-support threshold')
        ax.set_xticks([.8, .85, .9, .95, .99])
        ax.set_ylim(50, 101)
        ax.grid(axis='y', alpha=.2)
        ax.axhline(50, lw=.8, color='#888888', ls='--')
    axes[0].set_ylabel('Claims lacking resolved support (%)')
    axes[0].legend(frameon=False, loc='lower right', fontsize=8)
    fig.suptitle('Member-State short-window support across thresholds', fontweight='bold', fontsize=12)
    fig.supxlabel('Solid: stochastic-only lack of concordant resolution; dotted: final unresolved (curves overlap).', fontsize=8)
    for suffix in ('png', 'pdf'):
        fig.savefig(HERE / f'threshold_survival.{suffix}', bbox_inches='tight')
    plt.close(fig)
    claims = read(HERE / 'claim_level_continuous_support.csv')
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.5), sharey=True, layout='constrained')
    for ax, metric in zip(axes, METRICS):
        for scenario in SCALES:
            vals = np.sort([float(r['weakest_variant_support']) for r in claims if r['population'] == POPS[0] and r['scenario'] == scenario and r['metric'] == metric and r['weakest_variant_support'] != ''])
            ax.step(vals, np.arange(1, len(vals)+1)/len(vals), where='post', color=colors[scenario], lw=1.8,
                    label=f'{names[scenario]} (n={len(vals)})')
        ax.set_title(labels[metric], fontweight='bold')
        ax.set_xlabel('Weakest-variant sign support')
        ax.set_xlim(.4, 1.005)
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=.2)
        ax.legend(frameon=False, fontsize=7.5, loc='lower right')
    axes[0].set_ylabel('Cumulative fraction of concordant claims')
    fig.suptitle('Continuous support under paired amplitude scenarios', fontweight='bold', fontsize=12)
    fig.supxlabel('Conditional on concordant preferred classes; conflicts and missing coverage are excluded and reported separately.', fontsize=8)
    for suffix in ('png', 'pdf'):
        fig.savefig(HERE / f'continuous_support_distribution.{suffix}', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('baseline', 'run', 'aggregate', 'verify', 'report'))
    stage = parser.parse_args().stage
    {'baseline': baseline, 'run': run_scenarios, 'aggregate': aggregate, 'verify': verify_serialized, 'report': report}[stage]()
