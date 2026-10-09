"""Independent WP1B-R artifact checks and deterministic full rerun."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t

from wp1b_r import GEOS, MS, OUT, ROOT, VARIANTS, admissible, cmur, load_panel

FILES = ["branchA_medium_term.csv", "branchB1_anchor_medium_term.csv", "b2r_residual_pool.csv",
         "b2r_short_window_support.csv", "temporal_stability_comparison.csv"]

# Digests of claim-bearing and unaffected projections from the validated
# pre-correction outputs. These prevent a byte-identical rerun from masking
# an unintended change while allowing the bounded midpoint correction.
BASELINE = {
    "stochastic_rows": "2ff54632f9c7c312a8d47233a0decf637345c6687ec1af276933351eb9a6dda6",
    "turning_rows": "a2cc15f29445e113656c092612549f1107d61e4e1c5435ae19a4f2770cb84846",
    "deterministic_non_acceleration": "5a2f640e4f1c40c29e6e28af9946a52252bacde8b9d73ac0e95db9e550e9ecc0",
    "deterministic_noncenter_acceleration": "709375679d2ac4d333981320383a07f66c70f1ae185cf37e9e948ef0ef6928d9",
    "final_source_aware_classes": "6ff9a095b1d630131b9fa55cdb52360286b83384b2f1ca1f1cada22925637095",
    "final_stability_classes": "e055239a9a22d7574a40242cd89d67390a5203d7d749589dbf66cf4c2d0f8110",
    "medium_A": "37887c34fcc727f972880e918f85c24727e6b12b4c1c74b3270f6e0cacd27d6c",
    "medium_B1": "7092913b2e00236afee2b89c8633e847be333c3e2b107aaaacc645b62c174a8e",
    "residual_pool": "862cb601301dca5d95a98ebfb1f1c2fecacf45f3ab60cd3881ca2bd7ef67359e",
}


def projection_sha256(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def rows(name):
    with (OUT/name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def hashes():
    return {n:hashlib.sha256((OUT/n).read_bytes()).hexdigest() for n in FILES}


def main():
    failures, checks = [], {}
    def test(name, condition):
        checks[name] = bool(condition)
        if not condition: failures.append(name)
    p,a = load_panel()
    gate = json.loads((ROOT/"logs"/"wp1b_r_verification.json").read_text(encoding="utf-8"))
    A,B,R,S,C = (rows(n) for n in FILES)
    test("k3_gate_preceded_analysis",gate.get("gate")=="PASS" and gate.get("analysis_status")=="PASS" and all(gate["checks"].values()))
    test("medium_unique_exact_population",len(A)==len(B)==28 and {r["geo"] for r in A}==set(GEOS) and {r["geo"] for r in B}==set(GEOS))
    test("medium_fixed_2010_2022",all(r["start_year"]=="2010" and r["end_year"]=="2022" for r in A+B))
    test("branchA_all_annual_observations",all(int(r["n"])==13 for r in A))
    test("branchB1_verified_anchors_only",all(
        [int(y) for y in r["anchor_years"].split(";")] == [y for y in range(2010,2023) if admissible(a,r["geo"],y)]
        and int(r["n"])==len(r["anchor_years"].split(";")) for r in B))
    ols_ok=True
    for branch,records in (("A",A),("B1",B)):
        for r in records:
            years=np.asarray([y for y in range(2010,2023) if branch=="A" or admissible(a,r["geo"],y)],dtype=float)
            values=np.asarray([cmur(p,r["geo"],int(y)) for y in years])
            x=years-years.mean();slope=float(np.dot(x,values-values.mean())/np.dot(x,x))
            residual=values-(values.mean()+slope*x)
            se=float(np.sqrt(np.dot(residual,residual)/(len(years)-2)/np.dot(x,x)))
            q=float(student_t.ppf(.975,len(years)-2))
            ols_ok &= math.isclose(slope,float(r["slope_pp_per_year"]),abs_tol=1e-12)
            ols_ok &= math.isclose(slope-q*se,float(r["ci_low"]),abs_tol=1e-12)
            ols_ok &= math.isclose(slope+q*se,float(r["ci_high"]),abs_tol=1e-12)
    test("medium_ols_slope_and_95_interval_reproducible",ols_ok)
    test("residual_count_and_unique_keys",len(R)==135 and len({(r["geo"],r["interior_anchor_year"]) for r in R})==135)
    test("residuals_anchor_only_even_and_formula",all(
        all(admissible(a,r["geo"],int(r[k])) for k in ("left_anchor_year","interior_anchor_year","right_anchor_year"))
        and int(r["interior_anchor_year"])%2==0
        and math.isclose(float(r["signed_residual_pp"]),
            cmur(p,r["geo"],int(r["interior_anchor_year"]))-
            (cmur(p,r["geo"],int(r["left_anchor_year"]))+cmur(p,r["geo"],int(r["right_anchor_year"])))/2,
            abs_tol=1e-12) for r in R))
    count=Counter(r["geo"] for r in R)
    test("pooled_ms_count_130",sum(count[g] for g in MS)==130 and all(count[g]==int(gate["residual_count_by_geography"][g]) for g in GEOS))
    expected_short={(g,y,m,v) for g in GEOS for m,years in (("annual_change",range(2011,2025)),("turning",range(2012,2025)),("acceleration",range(2012,2025))) for y in years for v in VARIANTS}
    short_keys={(r["geo"],int(r["year"]),r["metric"],r["variant"]) for r in S}
    test("short_unique_full_grid",len(S)==len(expected_short)==len(short_keys) and short_keys==expected_short)
    test("endpoint_source_aware_unresolved",all(r["source_aware_estimate"]=="" and r["variant_class"] in ("UNRESOLVED","NOT_APPLICABLE") for r in S if int(r["year"])>=2023))
    test("country_variant_only_when_n_ge_4",all((r["variant_class"]!="NOT_APPLICABLE")==(count[r["geo"]]>=4) for r in S if r["variant"]=="country_signed"))
    test("all_simulated_estimates_bounded",all(0<=float(r["source_aware_estimate"])<=100 for r in S if r["source_aware_estimate"]!="" and r["metric"]=="turning"))
    test("bootstrap_sign_support_threshold",all(
        r["variant_class"] == ("POSITIVE" if float(r["p_positive"])>=.95 else "NEGATIVE" if float(r["p_negative"])>=.95 else "UNRESOLVED")
        for r in S if r["variant"] in ("pooled_signed","country_signed","pooled_absolute_random_sign") and r["metric"] in ("annual_change","acceleration") and r["source_aware_estimate"]!=""))
    # Build the midpoint-centred grid from the frozen anchor map, independent
    # of the production second-difference expression and its floating signs.
    centered_keys={(g,str(y),"acceleration","deterministic_midpoint") for g in GEOS
                   for y in range(2012,2023,2) if admissible(a,g,y-2) and admissible(a,g,y)}
    centered_rows=[r for r in S if (r["geo"],r["year"],r["metric"],r["variant"]) in centered_keys]
    test("midpoint_centered_structural_count_164",len(centered_keys)==len(centered_rows)==164)
    test("midpoint_centered_exact_zero_and_unresolved",all(
        r["source_aware_estimate"]=="0.0" and r["variant_class"]=="UNRESOLVED"
        for r in centered_rows))
    test("stochastic_rows_unchanged_from_validated_baseline",projection_sha256(
        [r for r in S if r["variant"]!="deterministic_midpoint"])==BASELINE["stochastic_rows"])
    test("turning_rows_unchanged_from_validated_baseline",projection_sha256(
        [r for r in S if r["metric"]=="turning"])==BASELINE["turning_rows"])
    test("deterministic_annual_and_turning_rows_unchanged",projection_sha256(
        [r for r in S if r["variant"]=="deterministic_midpoint" and r["metric"]!="acceleration"])
        ==BASELINE["deterministic_non_acceleration"])
    test("deterministic_noncenter_acceleration_unchanged",projection_sha256(
        [r for r in S if r["variant"]=="deterministic_midpoint" and r["metric"]=="acceleration"
         and (r["geo"],r["year"],r["metric"],r["variant"]) not in centered_keys])
        ==BASELINE["deterministic_noncenter_acceleration"])
    test("medium_and_residual_files_unchanged",hashes()["branchA_medium_term.csv"]==BASELINE["medium_A"]
         and hashes()["branchB1_anchor_medium_term.csv"]==BASELINE["medium_B1"]
         and hashes()["b2r_residual_pool.csv"]==BASELINE["residual_pool"])
    expected_cmp={(g,"2010-2022","medium_term_slope") for g in GEOS} | {(g,str(y),m) for g in GEOS for m,years in (("annual_change",range(2011,2025)),("turning",range(2012,2025)),("acceleration",range(2012,2025))) for y in years}
    cmp_keys={(r["geo"],r["year"],r["metric"]) for r in C}
    test("comparison_unique_full_grid",len(C)==len(expected_cmp)==len(cmp_keys) and cmp_keys==expected_cmp)
    test("comparison_endpoint_unresolved",all(r["source_aware_class"]=="UNRESOLVED" and r["stability_class"]=="UNRESOLVED" for r in C if r["year"]!="2010-2022" and int(r["year"])>=2023))
    short_by=defaultdict(dict)
    for r in S: short_by[(r["geo"],r["year"],r["metric"])][r["variant"]]=r
    consensus_ok=True
    for row in C:
        if row["metric"]=="medium_term_slope": continue
        variants=short_by[(row["geo"],row["year"],row["metric"])]
        classes=[x["variant_class"] for x in variants.values() if x["variant_class"]!="NOT_APPLICABLE"]
        agree=len(classes)>=3 and len(set(classes))==1
        final=classes[0] if agree and classes[0]!="UNRESOLVED" else "UNRESOLVED"
        consensus_ok &= row["source_aware_class"]==final and row["sensitivity_agreement"]==str(agree)
    test("predeclared_variant_consensus_rule",consensus_ok)
    test("final_source_aware_classes_unchanged",projection_sha256(
        [(r["geo"],r["year"],r["metric"],r["source_aware_class"]) for r in C])
        ==BASELINE["final_source_aware_classes"])
    test("final_stability_classes_unchanged",projection_sha256(
        [(r["geo"],r["year"],r["metric"],r["stability_class"]) for r in C])
        ==BASELINE["final_stability_classes"])
    expected_headlines={
        "medium_term_slope": Counter({"AGREEMENT":15,"UNRESOLVED":12}),
        "annual_change": Counter({"AGREEMENT":10,"UNRESOLVED":314}),
        "turning": Counter({"AGREEMENT":2,"UNRESOLVED":295}),
        "acceleration": Counter({"AGREEMENT":1,"UNRESOLVED":296}),
    }
    test("member_state_headline_counts_unchanged",all(
        Counter(r["stability_class"] for r in C if r["metric"]==metric and r["geo"] in MS
                and (metric=="medium_term_slope" or int(r["year"])<=2022))==counts
        for metric,counts in expected_headlines.items()))
    test("eu_aggregate_headline_counts_unchanged",all(
        Counter(r["stability_class"] for r in C if r["metric"]==metric and r["geo"]=="EU27_2020"
                and int(r["year"])<=2022)==Counter({"UNRESOLVED":n})
        for metric,n in (("annual_change",12),("turning",11),("acceleration",11)))
        and next(r for r in C if r["metric"]=="medium_term_slope" and r["geo"]=="EU27_2020")["stability_class"]=="AGREEMENT")
    test("ms_and_eu_separate_complete",sum(r["geo"] in MS for r in A)==27 and sum(r["geo"]=="EU27_2020" for r in A)==1)
    before=hashes()
    proc=subprocess.run([sys.executable,str(ROOT/"scripts"/"wp1b_r.py"),"--analyze"],cwd=ROOT,capture_output=True,text=True)
    after=hashes()
    test("deterministic_full_rerun",proc.returncode==0 and before==after)
    gate = json.loads((ROOT/"logs"/"wp1b_r_verification.json").read_text(encoding="utf-8"))
    gate.update(verification_status="PASS" if not failures else "FAIL",post_analysis_checks=checks,
                post_analysis_failures=failures,deterministic_rerun_sha256=after)
    (ROOT/"logs"/"wp1b_r_verification.json").write_text(json.dumps(gate,indent=2),encoding="utf-8")
    print(gate["verification_status"],len(checks),"post-analysis checks",failures)
    if failures: raise SystemExit(1)


if __name__=="__main__":
    main()
