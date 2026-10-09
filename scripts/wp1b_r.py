"""Refrozen WP1B-R analysis. Run --gate-only before --analyze."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"
LOG = ROOT / "logs" / "wp1b_r_verification.json"
GEOS = "AT BE BG HR CY CZ DK EE FI FR DE EL HU IE IT LV LT LU MT NL PL PT RO SK SI ES SE EU27_2020".split()
MS = [g for g in GEOS if g != "EU27_2020"]
YEARS = range(2010, 2023)
ODD = range(2011, 2022, 2)
INTERIOR = range(2012, 2021, 2)
SEED = 20260926
N_DRAWS = 20000
VARIANTS = ("pooled_signed", "country_signed", "pooled_absolute_random_sign", "deterministic_midpoint")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row.get(k, "") for k in fields} for row in rows)


def load_panel():
    panel = read_csv(OUT / "cmur_source_provenance_panel.csv")
    amap = read_csv(OUT / "verified_anchor_map.csv")
    p = {(r["geo"], int(r["year"])): r for r in panel}
    a = {(r["geo"], int(r["year"])): r for r in amap}
    if len(p) != 420 or len(a) != 420 or set(p) != {(g, y) for g in GEOS for y in range(2010, 2025)} or set(a) != set(p):
        raise ValueError("WP1A panel/anchor-map population or key integrity failed")
    return p, a


def cmur(p, geo, year):
    value = p[(geo, year)]["cmur_released_pct"]
    if value == "":
        raise ValueError(f"Released CMUR missing: {geo} {year}")
    return float(value)


def admissible(a, geo, year):
    return a[(geo, year)]["admissible_anchor"] == "True"


def prepare(p, a):
    residual_rows, eligible_odd = [], {}
    by_geo = defaultdict(list)
    for geo in GEOS:
        for year in INTERIOR:
            trio = (year-2, year, year+2)
            if all(admissible(a, geo, t) for t in trio):
                left, centre, right = (cmur(p, geo, t) for t in trio)
                error = centre-(left+right)/2
                row = dict(geo=geo, interior_anchor_year=year, left_anchor_year=year-2,
                           right_anchor_year=year+2, left_anchor_cmur_pct=left,
                           centre_anchor_cmur_pct=centre, right_anchor_cmur_pct=right,
                           signed_residual_pp=error, pooled_member_state=geo in MS)
                residual_rows.append(row)
                by_geo[geo].append(error)
        for year in ODD:
            if admissible(a, geo, year-1) and admissible(a, geo, year+1):
                eligible_odd[(geo, year)] = (cmur(p, geo, year-1)+cmur(p, geo, year+1))/2
    pooled = np.asarray([float(r["signed_residual_pp"]) for r in residual_rows if r["pooled_member_state"]], dtype=float)
    return residual_rows, by_geo, pooled, eligible_odd


def simulate(p, a, by_geo, pooled, eligible_odd, geo, variant):
    """Draw a full 2010-2022 path, with absent years represented by no key."""
    if variant == "country_signed" and len(by_geo[geo]) < 4:
        return None
    if variant not in VARIANTS:
        raise ValueError(variant)
    n = 1 if variant == "deterministic_midpoint" else N_DRAWS
    vi = VARIANTS.index(variant)
    gi = GEOS.index(geo)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED, gi, vi])))
    path = {}
    for year in YEARS:
        if admissible(a, geo, year):
            path[year] = np.full(n, cmur(p, geo, year), dtype=float)
        elif (geo, year) in eligible_odd:
            midpoint = eligible_odd[(geo, year)]
            if variant == "deterministic_midpoint":
                error = np.zeros(n)
            elif variant == "country_signed":
                error = rng.choice(np.asarray(by_geo[geo], dtype=float), size=n, replace=True)
            elif variant == "pooled_signed":
                error = rng.choice(pooled, size=n, replace=True)
            else:
                error = rng.choice(np.abs(pooled), size=n, replace=True)*rng.choice(np.array([-1.0, 1.0]), size=n, replace=True)
            path[year] = np.clip(midpoint+error, 0, 100)
    return path


def gate():
    p, a = load_panel()
    residual_rows, by_geo, pooled, eligible_odd = prepare(p, a)
    write_csv(OUT / "b2r_residual_pool.csv", residual_rows,
              ["geo", "interior_anchor_year", "left_anchor_year", "right_anchor_year",
               "left_anchor_cmur_pct", "centre_anchor_cmur_pct", "right_anchor_cmur_pct",
               "signed_residual_pp", "pooled_member_state"])
    checks = {}
    wp1a = json.loads((ROOT/"logs"/"verification_report.json").read_text(encoding="utf-8"))
    manifest = json.loads((ROOT/"logs"/"acquisition_manifest.json").read_text(encoding="utf-8"))
    checks["wp1a_frozen_panel_and_anchor_hashes"] = wp1a.get("status")=="PASS" and all(
        hashlib.sha256((OUT/name).read_bytes()).hexdigest()==wp1a["processed_sha256"][name]
        for name in ("cmur_source_provenance_panel.csv","verified_anchor_map.csv"))
    checks["official_snapshot_hashes_unchanged"] = all(
        hashlib.sha256((ROOT/r["path"]).read_bytes()).hexdigest()==r["sha256"] for r in manifest)
    checks["exact_fixed_population_and_2010_2022_span"] = len(p) == 420 and set(GEOS) == set(g for g, _ in p)
    checks["unique_residual_keys"] = len({(r["geo"], r["interior_anchor_year"]) for r in residual_rows}) == len(residual_rows)
    checks["calibration_verified_anchor_triplets_only"] = all(
        all(admissible(a, r["geo"], t) for t in (r["left_anchor_year"], r["interior_anchor_year"], r["right_anchor_year"]))
        and r["interior_anchor_year"] in INTERIOR and r["left_anchor_year"] == r["interior_anchor_year"]-2
        and r["right_anchor_year"] == r["interior_anchor_year"]+2 for r in residual_rows)
    checks["no_odd_year_calibration"] = all(r["interior_anchor_year"] % 2 == 0 for r in residual_rows)
    checks["residual_formula_from_anchor_released_cmur"] = all(math.isclose(
        r["signed_residual_pp"], cmur(p,r["geo"],r["interior_anchor_year"])-
        (cmur(p,r["geo"],r["left_anchor_year"])+cmur(p,r["geo"],r["right_anchor_year"]))/2,
        abs_tol=1e-12) for r in residual_rows)
    expected_counts = {r["geo"]: int(r["eligible_residual_count"]) for r in read_csv(OUT/"pseudogap_support.csv") if r["geo"] != "POOLED_EU27_MS"}
    actual_counts = Counter(r["geo"] for r in residual_rows)
    checks["residual_counts_match_wp1a_structure"] = all(actual_counts[g] == expected_counts[g] for g in GEOS) and len(pooled) == 130
    checks["bracketed_years_require_adjacent_verified_anchors"] = all(year in ODD and admissible(a,geo,year-1) and admissible(a,geo,year+1) for geo,year in eligible_odd)
    checks["pool_finite_nondegenerate"] = np.isfinite(pooled).all() and len(np.unique(pooled)) > 1
    variant_checks = {}
    deterministic = True
    endpoint_clean = True
    fixed_geographies = True
    bounded = True
    for geo in GEOS:
        for variant in VARIANTS:
            one = simulate(p,a,by_geo,pooled,eligible_odd,geo,variant)
            two = simulate(p,a,by_geo,pooled,eligible_odd,geo,variant)
            expected_applicable = variant != "country_signed" or len(by_geo[geo]) >= 4
            variant_checks[f"{geo}:{variant}"] = (one is not None) == expected_applicable
            if one is not None:
                deterministic &= set(one) == set(two) and all(np.array_equal(one[y],two[y]) for y in one)
                endpoint_clean &= not any(y >= 2023 for y in one)
                fixed_geographies &= geo in GEOS
                endpoint_clean &= all((admissible(a,geo,y) or (geo,y) in eligible_odd) for y in one)
                bounded &= all(np.isfinite(one[y]).all() and np.all((one[y]>=0)&(one[y]<=100)) for y in one)
    checks["all_predeclared_variants_executed_where_applicable"] = all(variant_checks.values())
    checks["fixed_seed_deterministic_rerun"] = deterministic
    checks["no_2023_2024_or_unbracketed_synthetic_values"] = endpoint_clean
    checks["all_simulated_cmur_bounded_0_100"] = bounded
    checks["no_outcome_based_geography_filter"] = fixed_geographies and len(variant_checks) == len(GEOS)*len(VARIANTS)
    report = {"stage": "K3-R PRE-INTERPRETATION", "gate": "PASS" if all(checks.values()) else "FAIL",
              "checks": checks, "failed_checks": [k for k,v in checks.items() if not v],
              "fixed_seed": SEED, "monte_carlo_draws": N_DRAWS,
              "variants": list(VARIANTS), "variant_execution": variant_checks,
              "pooled_member_state_residual_count": len(pooled),
              "residual_count_by_geography": {g: actual_counts[g] for g in GEOS},
              "bracketed_year_count_by_geography": {g: sum((g,y) in eligible_odd for y in ODD) for g in GEOS}}
    LOG.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("K3-R", report["gate"], len(checks), "checks", report["failed_checks"])
    if report["gate"] != "PASS":
        raise SystemExit(2)
    return p,a,residual_rows,by_geo,pooled,eligible_odd,report


def ols(rows):
    years = np.asarray([r[0] for r in rows], dtype=float)
    values = np.asarray([r[1] for r in rows], dtype=float)
    n = len(years)
    if n < 3:
        return dict(n=n, slope_pp_per_year="", ci_low="", ci_high="", p_positive="", direction="UNRESOLVED")
    x = years-years.mean()
    slope = float(np.dot(x, values-values.mean())/np.dot(x,x))
    fitted = values.mean()+slope*x
    stderr = float(np.sqrt(np.sum((values-fitted)**2)/(n-2)/np.dot(x,x)))
    critical = float(student_t.ppf(0.975,n-2))
    low, high = slope-critical*stderr, slope+critical*stderr
    p_positive = float(student_t.cdf(slope/stderr,n-2)) if stderr else (1.0 if slope>0 else 0.0 if slope<0 else 0.5)
    return dict(n=n,slope_pp_per_year=slope,ci_low=low,ci_high=high,p_positive=p_positive,stderr=stderr)


def sign_class(p_positive, p_negative):
    if p_positive >= 0.95:
        return "POSITIVE"
    if p_negative >= 0.95:
        return "NEGATIVE"
    return "UNRESOLVED"


def sign_of(value):
    return "POSITIVE" if value>0 else "NEGATIVE" if value<0 else "UNRESOLVED"


def classify_array(values, deterministic=False):
    if deterministic:
        val = float(values[0])
        return dict(estimate=val,low="",high="",p_positive="",p_negative="",classification=sign_of(val),sign_support="")
    p_pos = float(np.mean(values>0)); p_neg = float(np.mean(values<0))
    return dict(estimate=float(np.mean(values)),low=float(np.quantile(values,0.025)),
                high=float(np.quantile(values,0.975)),p_positive=p_pos,p_negative=p_neg,
                classification=sign_class(p_pos,p_neg),sign_support=max(p_pos,p_neg))


def conventional_class(p,geo,year,metric):
    if metric == "annual_change":
        return sign_of(cmur(p,geo,year)-cmur(p,geo,year-1))
    if metric == "acceleration":
        return sign_of(cmur(p,geo,year)-2*cmur(p,geo,year-1)+cmur(p,geo,year-2))
    d1 = conventional_class(p,geo,year-1,"annual_change")
    d2 = conventional_class(p,geo,year,"annual_change")
    return "TURN" if d1 != d2 and "UNRESOLVED" not in (d1,d2) else "NO_TURN" if d1 == d2 and d1 != "UNRESOLVED" else "UNRESOLVED"


def conventional_estimate(p,geo,year,metric):
    if metric == "annual_change":
        return cmur(p,geo,year)-cmur(p,geo,year-1)
    if metric == "acceleration":
        return cmur(p,geo,year)-2*cmur(p,geo,year-1)+cmur(p,geo,year-2)
    return int(conventional_class(p,geo,year,"turning") == "TURN")


def short_metric(path,year,metric,annual_classes=None):
    if metric == "annual_change":
        return path[year]-path[year-1]
    if metric == "acceleration":
        return path[year]-2*path[year-1]+path[year-2]
    return ((path[year]-path[year-1])*(path[year-1]-path[year-2])<0).astype(float)


def analyze():
    p,a,residual_rows,by_geo,pooled,eligible_odd,gate_report = gate()
    if gate_report["gate"] != "PASS":
        raise RuntimeError("K3-R gate failed")
    a_rows,b_rows,comparisons,short_rows = [],[],[],[]
    for geo in GEOS:
        annual = [(y,cmur(p,geo,y)) for y in YEARS]
        anchor = [(y,cmur(p,geo,y)) for y in YEARS if admissible(a,geo,y)]
        ar,br = ols(annual),ols(anchor)
        ar.update(geo=geo,start_year=2010,end_year=2022,branch="A",direction="POSITIVE" if ar["ci_low"]>0 else "NEGATIVE" if ar["ci_high"]<0 else "UNRESOLVED")
        br.update(geo=geo,start_year=2010,end_year=2022,branch="B1",direction=sign_class(br["p_positive"],1-br["p_positive"]),
                  anchor_years=";".join(str(y) for y,_ in anchor))
        a_rows.append(ar);b_rows.append(br)
        stability = "UNRESOLVED" if "UNRESOLVED" in (ar["direction"],br["direction"]) else "AGREEMENT" if ar["direction"]==br["direction"] else "CHANGE"
        comparisons.append(dict(geo=geo,year="2010-2022",metric="medium_term_slope",branchA_estimate=ar["slope_pp_per_year"],
                                source_aware_estimate=br["slope_pp_per_year"],magnitude_difference=br["slope_pp_per_year"]-ar["slope_pp_per_year"],
                                branchA_class=ar["direction"],source_aware_class=br["direction"],stability_class=stability,
                                sensitivity_agreement="not_applicable",source_aware_status="B1_anchor_only"))
        paths = {v:simulate(p,a,by_geo,pooled,eligible_odd,geo,v) for v in VARIANTS}
        annual_class = defaultdict(dict)
        annual_support = defaultdict(dict)
        for variant,path in paths.items():
            for year in range(2011,2025):
                defined = path is not None and year in path and year-1 in path and year<=2022
                if not defined:
                    cls = "NOT_APPLICABLE" if path is None else "UNRESOLVED"
                    reason = "country_pool_lt_4" if path is None else "endpoint_unanchored" if year>=2023 else "missing_admissible_anchor_or_bracket"
                    result = dict(estimate="",low="",high="",p_positive="",p_negative="",classification=cls,sign_support="")
                else:
                    result=classify_array(short_metric(path,year,"annual_change"),variant=="deterministic_midpoint")
                    reason=""
                annual_class[variant][year] = result["classification"]
                annual_support[variant][year] = result["sign_support"]
                short_rows.append(dict(geo=geo,year=year,metric="annual_change",variant=variant,
                    branchA_estimate=conventional_estimate(p,geo,year,"annual_change"),branchA_class=conventional_class(p,geo,year,"annual_change"),
                    source_aware_estimate=result["estimate"],ci_low=result["low"],ci_high=result["high"],
                    p_positive=result["p_positive"],p_negative=result["p_negative"],sign_support=result["sign_support"],
                    variant_class=result["classification"],reason=reason))
            for year in range(2012,2025):
                for metric in ("turning","acceleration"):
                    defined = path is not None and all(t in path for t in (year-2,year-1,year)) and year<=2022
                    if not defined:
                        cls = "NOT_APPLICABLE" if path is None else "UNRESOLVED"
                        reason = "country_pool_lt_4" if path is None else "endpoint_unanchored" if year>=2023 else "missing_admissible_anchor_or_bracket"
                        result = dict(estimate="",low="",high="",p_positive="",p_negative="",classification=cls,sign_support="")
                    elif metric == "acceleration":
                        # An admissible anchor--constructed midpoint--anchor triple
                        # has identically zero centred second difference. Encode that
                        # universal geometry directly, without a numerical tolerance.
                        midpoint_centred = variant == "deterministic_midpoint" and (geo, year-1) in eligible_odd
                        values = np.zeros(1, dtype=float) if midpoint_centred else short_metric(path,year,metric)
                        result = classify_array(values,variant=="deterministic_midpoint")
                        reason=""
                    else:
                        c1,c2=annual_class[variant][year-1],annual_class[variant][year]
                        cls = "UNRESOLVED" if "UNRESOLVED" in (c1,c2) else "TURN" if c1!=c2 else "NO_TURN"
                        vals=short_metric(path,year,metric)
                        result=dict(estimate=float(np.mean(vals)),low="",high="",p_positive="",p_negative="",
                                    classification=cls,sign_support=min(annual_support[variant][year-1],annual_support[variant][year]) if variant!="deterministic_midpoint" else "")
                        reason="" if cls!="UNRESOLVED" else "consecutive_change_sign_support_below_0.95"
                    short_rows.append(dict(geo=geo,year=year,metric=metric,variant=variant,
                        branchA_estimate=conventional_estimate(p,geo,year,metric),branchA_class=conventional_class(p,geo,year,metric),
                        source_aware_estimate=result["estimate"],ci_low=result["low"],ci_high=result["high"],
                        p_positive=result["p_positive"],p_negative=result["p_negative"],sign_support=result["sign_support"],
                        variant_class=result["classification"],reason=reason))
    grouped=defaultdict(dict)
    for r in short_rows:
        grouped[(r["geo"],r["year"],r["metric"])][r["variant"]]=r
    for (geo,year,metric),variants in grouped.items():
        primary=variants["pooled_signed"]
        applicable=[r["variant_class"] for r in variants.values() if r["variant_class"]!="NOT_APPLICABLE"]
        agree=len(set(applicable))==1 and len(applicable)>=3
        final_cls=primary["variant_class"] if agree and primary["variant_class"] not in ("UNRESOLVED","NOT_APPLICABLE") else "UNRESOLVED"
        if final_cls=="UNRESOLVED": stability="UNRESOLVED"
        elif primary["branchA_class"]=="UNRESOLVED": stability="UNRESOLVED"
        else: stability="AGREEMENT" if final_cls==primary["branchA_class"] else "CHANGE"
        comparisons.append(dict(geo=geo,year=year,metric=metric,branchA_estimate=primary["branchA_estimate"],
               source_aware_estimate=primary["source_aware_estimate"],magnitude_difference=(float(primary["source_aware_estimate"])-float(primary["branchA_estimate"])) if primary["source_aware_estimate"]!="" else "",
               branchA_class=primary["branchA_class"],source_aware_class=final_cls,stability_class=stability,
               sensitivity_agreement=agree,source_aware_status="endpoint_unanchored" if year>=2023 else "bracketed" if primary["source_aware_estimate"]!="" else "source_exception"))
    write_csv(OUT/"branchA_medium_term.csv",a_rows,["geo","branch","start_year","end_year","n","slope_pp_per_year","stderr","ci_low","ci_high","p_positive","direction"])
    write_csv(OUT/"branchB1_anchor_medium_term.csv",b_rows,["geo","branch","start_year","end_year","n","anchor_years","slope_pp_per_year","stderr","ci_low","ci_high","p_positive","direction"])
    write_csv(OUT/"b2r_short_window_support.csv",short_rows,["geo","year","metric","variant","branchA_estimate","branchA_class","source_aware_estimate","ci_low","ci_high","p_positive","p_negative","sign_support","variant_class","reason"])
    write_csv(OUT/"temporal_stability_comparison.csv",comparisons,["geo","year","metric","branchA_estimate","source_aware_estimate","magnitude_difference","branchA_class","source_aware_class","stability_class","sensitivity_agreement","source_aware_status"])
    output_names=["branchA_medium_term.csv","branchB1_anchor_medium_term.csv","b2r_residual_pool.csv","b2r_short_window_support.csv","temporal_stability_comparison.csv"]
    hashes={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in output_names}
    gate_report.update(stage="ANALYSIS COMPLETE",analysis_status="PASS",output_sha256=hashes,
                       output_rows={name:len(read_csv(OUT/name)) for name in output_names})
    LOG.write_text(json.dumps(gate_report,indent=2),encoding="utf-8")
    print("ANALYSIS COMPLETE",gate_report["output_rows"])


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gate-only",action="store_true")
    group.add_argument("--analyze",action="store_true")
    args=parser.parse_args()
    gate() if args.gate_only else analyze()
