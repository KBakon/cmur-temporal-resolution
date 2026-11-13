"""Independent serialized-cell checks and deterministic supplement builds only."""
from __future__ import annotations
from pathlib import Path
from collections import Counter
from statistics import median
import csv, hashlib, json, math, re, subprocess, sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
MAN = ROOT / "manuscript"
NAMES = ("tableS1_country_slopes.tex", "tableS2_country_intervals.tex",
         "tableS3_threshold_grid.tex", "tableS4_continuous_support.tex", "tableS5_diagnostics.tex")
OUTPUTS = tuple(MAN / "supplementary" / n for n in NAMES)
RULES = ("COMMON_TWO_SIDED_95", "COMMON_DIRECTIONAL_95")
SCENARIOS = ("UNSCALED_STRESS", "BROWNIAN_SD_1_SQRT2")
POPS = ("Member States", "EU27_2020")
METRICS = ("annual_change", "turning", "acceleration")
GRID = (.80, .85, .90, .925, .95, .975, .99)
LABELS = dict(annual_change="Annual change", turning="Turning", acceleration="Acceleration")

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def inventory():
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob("*") if p.is_file() and p not in OUTPUTS}

def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))

def cells(path, prefix=None):
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if " & " in line and line.endswith(r" \\") and (prefix is None or line.startswith(prefix)):
            result.append([v.strip().replace(r"\_", "_") for v in line[:-3].split("&")])
    return result

def main():
    before = inventory()
    runs = []
    for _ in range(2):
        proc = subprocess.run([sys.executable, "-B", str(MAN/"scripts/build_supplement_assets.py")],
                              cwd=ROOT, capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(proc.stdout+proc.stderr)
        runs.append(json.loads(proc.stdout))
    checks = {}
    def check(name, ok):
        checks[name] = bool(ok)
        if not ok:
            raise AssertionError(name)
    check("two_builds_byte_identical", runs[0]["output_sha256"] == runs[1]["output_sha256"])
    check("all_other_repository_paths_unchanged", before == inventory())
    freeze = (ROOT/"artifacts/SCIENTIFIC_FREEZE_HANDOFF.md").read_text(encoding="utf-8")
    source = {}
    for i in range(1, 14):
        if i == 8:
            continue
        key = f"N{i:02d}"
        line = next(l for l in freeze.splitlines() if l.startswith("| "+key+" |"))
        cols = line.split("|")
        path = ROOT/cols[2].strip().strip(chr(96))
        check(key+"_hash", digest(path) == cols[-2].strip().strip(chr(96)))
        source[key] = json.loads(path.read_text()) if path.suffix == ".json" else rows(path)
    primary = {r["geo"]: r for r in source["N01"] if r["rule"] == RULES[0]}
    secondary = {r["geo"]: r for r in source["N01"] if r["rule"] == RULES[1]}
    geos = set(primary)
    def country_rows(p):
        return [r for r in cells(p) if r[0] in geos]
    s1 = country_rows(OUTPUTS[0])
    check("S1_28_unique_geographies", len(s1) == len({r[0] for r in s1}) == 28 and {r[0] for r in s1} == geos)
    for row in s1:
        p, s = primary[row[0]], secondary[row[0]]
        expected = [p["geo"], f"{float(p['A_slope_pp_per_year']):.4f}", f"{float(p['B1_slope_pp_per_year']):.4f}", p["B1_n"],
                    p["A_class"], p["B1_class"], p["stability_class"], s["A_class"], s["B1_class"], s["stability_class"]]
        check("S1_cells_"+row[0], row == expected)
    for index, name, lock in ((6, "primary", (12, 0, 15)), (9, "secondary", (15, 0, 12))):
        c = Counter(r[index] for r in s1 if r[0] != "EU27_2020")
        check("S1_"+name+"_MS", tuple(c[k] for k in ("AGREEMENT", "CHANGE", "UNRESOLVED")) == lock)
        check("S1_"+name+"_EU", next(r for r in s1 if r[0] == "EU27_2020")[index] == "AGREEMENT")
    s2 = country_rows(OUTPUTS[1])
    check("S2_28_unique_geographies", len(s2) == len({r[0] for r in s2}) == 28)
    wa, wb = [], []
    for row in s2:
        p = primary[row[0]]
        a = float(p["A_ci_high"])-float(p["A_ci_low"])
        b = float(p["B1_ci_high"])-float(p["B1_ci_low"])
        expected = [p["geo"], f"[{float(p['A_ci_low']):.4f}, {float(p['A_ci_high']):.4f}]",
                    f"[{float(p['B1_ci_low']):.4f}, {float(p['B1_ci_high']):.4f}]", f"{a:.4f}", f"{b:.4f}",
                    f"{b/a:.4f}" if a else "--", f"{float(p['A_p_positive']):.6f}", f"{float(p['B1_p_positive']):.6f}"]
        check("S2_cells_"+row[0], row == expected)
        if row[0] != "EU27_2020":
            wa.append(a); wb.append(b)
    check("S2_27_wider", sum(b>a for a,b in zip(wa, wb)) == 27)
    check("S2_medians", (f"{median(wa):.4f}", f"{median(wb):.4f}", f"{median(b/a for a,b in zip(wa,wb)):.3f}") == ("0.3386","0.6219","1.858"))
    s2text = OUTPUTS[1].read_text()
    check("S2_summary_note", all(s in s2text for s in ("0.3386", "0.6219", "1.858", "27/27")))
    s3 = cells(OUTPUTS[2], ("Member States &", r"EU27\_2020 &"))
    key3 = lambda r: (r[0], r[1], r[2], float(r[3]))
    check("S3_84_unique_groups", len(s3) == len({key3(r) for r in s3}) == 84)
    expected_grid = {(p, LABELS[m], s, t) for p in POPS for m in METRICS for s in ("U","B") for t in GRID}
    check("S3_exact_grid", {key3(r) for r in s3} == expected_grid)
    for row in s3:
        sc = SCENARIOS[("U","B").index(row[2])]
        metric = next(m for m, label in LABELS.items() if label == row[1])
        n06 = [r for r in source["N06"] if (r["population"],r["metric"],r["scenario"],float(r["threshold"])) == (row[0],metric,sc,float(row[3]))]
        check("S3_source_unique_"+str(key3(row)),len(n06)==1)
        r = n06[0]
        values = [r[k] for k in ("agreement","resolved_class_disagreement","unresolved","total")]
        check("S3_cells_"+str(key3(row)), row[4:8] == values and row[8] == f"{100*int(r['unresolved'])/int(r['total']):.2f}")
        claims = [r for r in source["N07"] if (r["population"],r["metric"],r["scenario"],float(r["threshold"])) == (row[0],metric,sc,float(row[3]))]
        c = Counter(r["stability_class"] for r in claims)
        check("S3_claim_recount_"+str(key3(row)), tuple(c[k] for k in ("AGREEMENT","CHANGE","UNRESOLVED")) == tuple(map(int,row[4:7])) and len(claims)==int(row[7]))
        if row[0] == POPS[0]:
            check("S3_majority_"+str(key3(row)), int(row[6])/int(row[7]) > .5)
    be = [r for r in source["N07"] if (r["geo"],r["year"],r["metric"],r["scenario"]) == ("BE","2019","annual_change",SCENARIOS[1]) and float(r["threshold"]) in (.8,.85)]
    check("BE_2019_B_adverse_case", len(be)==2 and all((r["branchA_class"],r["source_aware_class"],r["stability_class"]) == ("NEGATIVE","POSITIVE","CHANGE") for r in be))
    s4 = cells(OUTPUTS[3], ("Member States &", r"EU27\_2020 &"))
    check("S4_12_unique_groups", len(s4)==len({tuple(r[:3]) for r in s4})==12)
    # NumPy is used solely to verify frozen descriptive quantiles from stored claim support.
    import numpy as np
    for row in s4:
        sc=SCENARIOS[("U","B").index(row[2])]
        m=next(m for m,label in LABELS.items() if label==row[1])
        summaries=[r for r in source["N10"] if (r["population"],r["metric"],r["scenario"],r["subset"]) == (row[0],m,sc,"scenario_concordant")]
        check("S4_summary_unique_"+str(row[:3]),len(summaries)==1)
        r=summaries[0]
        check("S4_cells_"+str(row[:3]), row[3:8]==[r[k] for k in ("total_claims","concordant_with_support","conflict","preferred_tie","missing_coverage")] and row[8:]==[f"{float(r[k]):.3f}" for k in ("p10","q25","median","q75","p90")])
        claims=[x for x in source["N09"] if (x["population"],x["metric"],x["scenario"]) == (row[0],m,sc)]
        c=Counter(x["preferred_class_status"] for x in claims)
        values=[float(x["weakest_variant_support"]) for x in claims if x["preferred_class_status"]=="CONCORDANT"]
        check("S4_denominators_"+str(row[:3]), list(map(int,row[3:8]))==[len(claims),len(values),c["CONFLICT"],c["PREFERRED_TIE"],c["MISSING_COVERAGE"]])
        check("S4_quantiles_N09_"+str(row[:3]),row[8:]==[f"{x:.3f}" for x in np.quantile(values,[.1,.25,.5,.75,.9])])
    check("S4_caption_semantics",all(x in OUTPUTS[3].read_text() for x in ("conditional, descriptive","need not be paired","not posterior, joint or calibrated")))
    s5 = cells(OUTPUTS[4])
    mid=[r for r in s5 if r[0].startswith(("Literal omission","Exploratory stochastic-only"))]
    check("S5_midpoint_four_rows",len(mid)==4 and all(r[2:]==["0","0"] for r in mid))
    for row in mid:
        diag="LOO_deterministic_midpoint" if row[0].startswith("Literal") else "EXPLORATORY_midpoint_structural_reference"
        rr=[r for r in source["N04"] if (r["diagnostic"],r["population"],r["metric"])==(diag,row[1],"acceleration") and int(r["year"])<=2022]
        check("S5_midpoint_counts_"+str(row[:2]), len(rr)==(297 if row[1]==POPS[0] else 11) and sum(r["newly_source_aware_resolved"]=="True" for r in rr)==int(row[2]) and sum(r["newly_comparison_resolved"]=="True" for r in rr)==int(row[3]))
    pairs=[r for r in source["N12"] if r["population"]==POPS[0]]
    c=Counter(r["ordering"] for r in pairs)
    check("S5_pairing_76_26_54_20_2",(len(pairs),len({r["geo"] for r in pairs}),c["8_GREATER"],c["4_GREATER"],c["EQUAL"])==(76,26,54,20,2))
    expected={"Pairs":"76","Countries":"26","Median absolute four-year residual (pp)":"0.525","Median absolute eight-year residual (pp)":"0.750","Eight-year absolute residual larger":"54","Four-year absolute residual larger":"20","Equal absolute residual":"2"}
    for label,value in expected.items():
        rr=[r for r in s5 if r[0]==label]
        check("S5_display_"+label,len(rr)==1 and rr[0][1]==value)
    for span,value in (("4",.525),("8",.750)):
        rr=[r for r in source["N13"] if (r["population"],r["paired_same_geo_centre"],r["outer_span_years"])==(POPS[0],"True",span)]
        summary=[r for r in source["N11"] if (r["population"],r["subset"],r["outer_span_years"])==(POPS[0],"paired_same_geo_centre",span)]
        check("S5_N11_N13_"+span,len(summary)==1 and len(rr)==76 and round(median(float(r["absolute_residual"]) for r in rr),3)==round(float(summary[0]["median_absolute_residual"]),3)==value)
    maintex=(MAN/"main.tex").read_text()
    check("all_five_inputs_wired", all("supplementary/"+name[:-4] in maintex for name in NAMES))
    check("no_other_files_changed_final",before==inventory())
    print(json.dumps(dict(all_supplement_checks_pass=True,check_count=len(checks),checks=checks,
                         source_sha256=runs[-1]["source_sha256"],output_sha256=runs[-1]["output_sha256"],
                         two_builds_byte_identical=True,scientific_analysis_rerun=False),indent=2))

if __name__=="__main__":
    main()
