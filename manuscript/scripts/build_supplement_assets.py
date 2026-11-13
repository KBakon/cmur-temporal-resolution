"""Supplement Tables S1-S5 from frozen serialized Number Authority only."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
MAN = ROOT / "manuscript"
SUP = MAN / "supplementary"
NAMES = (
    "tableS1_country_slopes.tex", "tableS2_country_intervals.tex",
    "tableS3_threshold_grid.tex", "tableS4_continuous_support.tex",
    "tableS5_diagnostics.tex",
)
OUTPUTS = tuple(SUP / name for name in NAMES)
RULES = ("COMMON_TWO_SIDED_95", "COMMON_DIRECTIONAL_95")
SCENARIOS = ("UNSCALED_STRESS", "BROWNIAN_SD_1_SQRT2")
POPS = ("Member States", "EU27_2020")
METRICS = ("annual_change", "turning", "acceleration")
GRID = (.80, .85, .90, .925, .95, .975, .99)
LABELS = dict(annual_change="Annual change", turning="Turning", acceleration="Acceleration")
END = r" \\"

def require(ok, message):
    if not ok:
        raise ValueError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))

def source_data():
    freeze = (ROOT / "artifacts/SCIENTIFIC_FREEZE_HANDOFF.md").read_text(encoding="utf-8")
    paths, hashes, data = {}, {}, {}
    for i in range(1, 14):
        if i == 8:  # U/B change rows are not required by these supplementary tables.
            continue
        key = f"N{i:02d}"
        lines = [line for line in freeze.splitlines() if line.startswith("| " + key + " |")]
        require(len(lines) == 1, f"Missing/duplicate freeze index {key}")
        cells = lines[0].split("|")
        relative = cells[2].strip().strip(chr(96))
        sha = cells[-2].strip().strip(chr(96))
        path = ROOT / relative
        require(path.is_file() and digest(path) == sha, f"SCIENTIFIC_AUTHORITY_CONFLICT: {key} hash")
        paths[key], hashes[key] = path, sha
        data[key] = json.loads(path.read_text(encoding="utf-8")) if path.suffix == ".json" else csv_rows(path)
    return data, paths, hashes

def unique(rows, key):
    result = {}
    for row in rows:
        k = key(row)
        require(k not in result, f"Duplicate group {k}")
        result[k] = row
    return result

def adu(rows):
    c = Counter(r["stability_class"] for r in rows)
    return tuple(c[k] for k in ("AGREEMENT", "CHANGE", "UNRESOLVED"))

def tex(value):
    return str(value).replace("_", r"\_")

def line(values):
    return " & ".join(map(tex, values)) + END

def longtable(label, caption, headers, specs, rows, notes="", size=r"\footnotesize", padding=3):
    # Headers are trusted LaTeX; escaping underscores would corrupt math subscripts.
    header = " & ".join(headers) + END
    return "\n".join([
        r"\begingroup", size, rf"\setlength{{\tabcolsep}}{{{padding}pt}}",
        r"\setlength{\LTcapwidth}{\linewidth}",
        r"\renewcommand{\arraystretch}{1.1}", rf"\begin{{longtable}}{{{specs}}}",
        rf"\caption{{{caption}}}\label{{tab:{label}}}" + END,
        r"\toprule", header, r"\midrule", r"\endfirsthead",
        rf"\multicolumn{{{len(headers)}}}{{l}}{{\textit{{Table {label} continued}}}}" + END,
        r"\toprule", header, r"\midrule", r"\endhead",
        r"\midrule", rf"\multicolumn{{{len(headers)}}}{{r}}{{\textit{{Continued on next page}}}}" + END,
        r"\endfoot", r"\bottomrule", r"\endlastfoot",
        *[line(row) for row in rows], r"\end{longtable}",
        notes, r"\endgroup", "",
    ])

def prepare():
    data, paths, hashes = source_data()
    n01, n05, n06, n07 = [data[k] for k in ("N01", "N05", "N06", "N07")]
    med = unique(n01, lambda r: (r["rule"], r["geo"]))
    require(len(med) == 56, "N01 complete 56-row grid")
    geos = sorted({r["geo"] for r in n01 if r["population"] == POPS[0]})
    require(len(geos) == 27 and "EU27_2020" not in geos, "Member-State population")
    geos += ["EU27_2020"]
    for rule, target in zip(RULES, ((12, 0, 15), (15, 0, 12))):
        require(set(g for ru, g in med if ru == rule) == set(geos), "N01 missing geography")
        for pop in POPS:
            rr = [r for r in n01 if r["rule"] == rule and r["population"] == pop]
            wanted = target if pop == POPS[0] else (1, 0, 0)
            require(adu(rr) == wanted, "Medium-term count lock")
            check = [r for r in n05 if r["rule"] == rule and r["population"] == pop]
            require(len(check) == 1, "N05 group missing/duplicate")
            require(tuple(int(check[0][k]) for k in ("agreement", "resolved_class_disagreement", "unresolved")) == wanted, "N05 mismatch")
    rows1, rows2, wa, wb = [], [], [], []
    for geo in geos:
        a, b = [med[rule, geo] for rule in RULES]
        require(all(a[k] == b[k] for k in a if k.startswith(("A_", "B1_")) and not k.endswith(("class", "changed"))), "Fixed fits differ")
        rows1.append([geo, f"{float(a['A_slope_pp_per_year']):.4f}", f"{float(a['B1_slope_pp_per_year']):.4f}", a["B1_n"],
                      a["A_class"], a["B1_class"], a["stability_class"], b["A_class"], b["B1_class"], b["stability_class"]])
        x = float(a["A_ci_high"]) - float(a["A_ci_low"])
        y = float(a["B1_ci_high"]) - float(a["B1_ci_low"])
        ratio = f"{y/x:.4f}" if x != 0 else "--"
        rows2.append([geo, f"[{float(a['A_ci_low']):.4f}, {float(a['A_ci_high']):.4f}]",
                      f"[{float(a['B1_ci_low']):.4f}, {float(a['B1_ci_high']):.4f}]",
                      f"{x:.4f}", f"{y:.4f}", ratio, f"{float(a['A_p_positive']):.6f}", f"{float(a['B1_p_positive']):.6f}"])
        if geo != "EU27_2020":
            wa.append(x); wb.append(y)
    precision = (median(wa), median(wb), median(y/x for x, y in zip(wa, wb)), sum(y > x for x, y in zip(wa, wb)))
    require(tuple(round(v, n) for v, n in zip(precision[:3], (4, 4, 3))) == (.3386, .6219, 1.858) and precision[3] == 27, "Precision lock")
    for r in n05:
        if r["population"] == POPS[0]:
            for actual, name in zip(precision[:3], ("median_A_CI_width", "median_B1_CI_width", "median_paired_B1_A_CI_width_ratio")):
                require(math.isclose(actual, float(r[name]), rel_tol=0, abs_tol=1e-12), "N05 precision mismatch")

    threshold = unique(n06, lambda r: (r["population"], r["metric"], r["scenario"], float(r["threshold"])))
    expected = {(p, m, s, t) for p in POPS for m in METRICS for s in SCENARIOS for t in GRID}
    require(set(threshold) == expected, "N06 complete threshold grid")
    detail = {}
    require(len(n07) == 13328, "N07 row count")
    for r in n07:
        detail.setdefault((r["population"], r["metric"], r["scenario"], float(r["threshold"])), []).append(r)
    require(set(detail) == expected, "N07 complete grid")
    rows3 = []
    locks = {SCENARIOS[0]: ((10, 0, 314), (2, 0, 295), (1, 0, 296)),
             SCENARIOS[1]: ((29, 0, 295), (13, 0, 284), (7, 0, 290))}
    for pop in POPS:
        for metric in METRICS:
            for scenario, name in zip(SCENARIOS, ("U", "B")):
                for tau in GRID:
                    key = (pop, metric, scenario, tau)
                    r, rr = threshold[key], detail[key]
                    total = (324 if metric == METRICS[0] else 297) if pop == POPS[0] else (12 if metric == METRICS[0] else 11)
                    counts = tuple(int(r[k]) for k in ("agreement", "resolved_class_disagreement", "unresolved"))
                    require(len(rr) == total == int(r["total"]) and len(unique(rr, lambda x: (x["geo"], x["year"], x["metric"]))) == total, "N07 key integrity")
                    require(adu(rr) == counts and sum(counts) == total, "N06/N07 mismatch")
                    require(math.isclose(counts[2]/total, float(r["unresolved_proportion"]), abs_tol=1e-15), "N06 share")
                    if pop == POPS[0]:
                        require(counts[2] > total/2, "Majority unresolved lock")
                    if tau == .95:
                        require(counts == (locks[scenario][METRICS.index(metric)] if pop == POPS[0] else (0, 0, total)), "Reference threshold lock")
                    rows3.append([pop, LABELS[metric], name, f"{tau:.3f}".rstrip("0") if tau in (.925, .975) else f"{tau:.2f}",
                                  *counts, total, f"{100*counts[2]/total:.2f}"])
    adverse = [r for r in n07 if r["geo"] == "BE" and r["year"] == "2019" and r["metric"] == METRICS[0] and r["scenario"] == SCENARIOS[1] and float(r["threshold"]) in (.8, .85)]
    require(len(adverse) == 2 and all((r["branchA_class"], r["source_aware_class"], r["stability_class"]) == ("NEGATIVE", "POSITIVE", "CHANGE") for r in adverse), "BE adverse lock")

    summaries = unique(data["N10"], lambda r: (r["scenario"], r["population"], r["metric"], r["subset"]))
    claims = unique(data["N09"], lambda r: (r["scenario"], r["geo"], r["year"], r["metric"]))
    require(len(claims) == 1904, "N09 complete claim grid")
    rows4 = []
    for pop in POPS:
        for metric in METRICS:
            for scenario, name in zip(SCENARIOS, ("U", "B")):
                r = summaries.get((scenario, pop, metric, "scenario_concordant"))
                require(r is not None, "Missing scenario-concordant summary")
                rr = [x for x in claims.values() if (x["scenario"], x["population"], x["metric"]) == (scenario, pop, metric)]
                states = Counter(x["preferred_class_status"] for x in rr)
                vals = sorted(float(x["weakest_variant_support"]) for x in rr if x["preferred_class_status"] == "CONCORDANT")
                require(len(rr) == int(r["total_claims"]) and len(vals) == int(r["concordant_with_support"]), "Continuous denominator")
                for field, status in (("conflict", "CONFLICT"), ("preferred_tie", "PREFERRED_TIE"), ("missing_coverage", "MISSING_COVERAGE")):
                    require(int(r[field]) == states[status], "Continuous subset mismatch")
                require(math.isclose(median(vals), float(r["median"]), abs_tol=1e-15), "Continuous median mismatch")
                require(sum(states.values()) == int(r["total_claims"]), "Subset accounting")
                rows4.append([pop, LABELS[metric], name, r["total_claims"], r["concordant_with_support"], r["conflict"], r["preferred_tie"], r["missing_coverage"],
                              *[f"{float(r[k]):.3f}" for k in ("p10", "q25", "median", "q75", "p90")]])
    for scenario, medians in zip(SCENARIOS, ((.632, .602, .557), (.709, .649, .567))):
        require(tuple(round(float(summaries[scenario, POPS[0], m, "scenario_concordant"]["median"]), 3) for m in METRICS) == medians, "C7 median locks")

    rows5a = []
    for diag, name in (("LOO_deterministic_midpoint", "Literal omission (minimum three)"),
                       ("EXPLORATORY_midpoint_structural_reference", "Exploratory stochastic-only (minimum two)")):
        for pop in POPS:
            rr = [r for r in data["N04"] if r["diagnostic"] == diag and r["population"] == pop and r["metric"] == "acceleration" and int(r["year"]) <= 2022]
            require(len(rr) == (297 if pop == POPS[0] else 11), "Midpoint complete acceleration grid")
            new_source = sum(r["newly_source_aware_resolved"] == "True" for r in rr)
            new_final = sum(r["newly_comparison_resolved"] == "True" for r in rr)
            require(new_source == new_final == 0, "Midpoint omission lock")
            summary = [r for r in data["N02"]["influence"] if r["diagnostic"] == diag and r["population"] == pop and r["metric"] == "acceleration" and r["period"] == "2010-2022"]
            require(len(summary) == 1 and summary[0]["newly_source_aware_resolved"] == new_source and summary[0]["newly_comparison_resolved"] == new_final, "N02 influence mismatch")
            rows5a.append([name, pop, new_source, new_final])
    for pop in POPS:
        rr = [r for r in data["N03"] if r["population"] == pop and r["metric"] == "acceleration" and int(r["year"]) <= 2022]
        summary = next(r for r in data["N02"]["causes"] if r["population"] == pop and r["metric"] == "acceleration" and r["period"] == "2010-2022")
        require(len(rr) == summary["unresolved"] and sum(r["midpoint_only_veto"] == "True" for r in rr) == 0, "Midpoint cause reconciliation")
    pairrows = [r for r in data["N12"] if r["population"] == POPS[0]]
    require(len(unique(pairrows, lambda r: (r["geo"], r["centre_year"]))) == 76 and len({r["geo"] for r in pairrows}) == 26, "76/26 pairing")
    order = Counter(r["ordering"] for r in pairrows)
    require(tuple(order[k] for k in ("8_GREATER", "4_GREATER", "EQUAL")) == (54, 20, 2), "Pair ordering")
    spans = []
    for span, expected in (("4", .525), ("8", .750)):
        candidates = [r for r in data["N11"] if (r["population"], r["subset"], r["outer_span_years"]) == (POPS[0], "paired_same_geo_centre", span)]
        require(len(candidates) == 1, "N11 pair summary missing/duplicate")
        r = candidates[0]
        cases = [x for x in data["N13"] if x["population"] == POPS[0] and x["paired_same_geo_centre"] == "True" and x["outer_span_years"] == span]
        require(len(unique(cases, lambda x: (x["geo"], x["centre_year"]))) == int(r["eligible_count"]) == 76, "N13 paired cases")
        require({(x["geo"], x["centre_year"]) for x in cases} == {(x["geo"], x["centre_year"]) for x in pairrows}, "Pair keys mismatch")
        require(round(median(float(x["absolute_residual"]) for x in cases), 3) == round(float(r["median_absolute_residual"]), 3) == expected, "Span median")
        spans.append(r)
    rows5b = [["Pairs", len(pairrows)], ["Countries", len({r["geo"] for r in pairrows})],
              ["Median absolute four-year residual (pp)", f"{float(spans[0]['median_absolute_residual']):.3f}"],
              ["Median absolute eight-year residual (pp)", f"{float(spans[1]['median_absolute_residual']):.3f}"],
              ["Eight-year absolute residual larger", order["8_GREATER"]],
              ["Four-year absolute residual larger", order["4_GREATER"]], ["Equal absolute residual", order["EQUAL"]]]
    return dict(rows1=rows1, rows2=rows2, rows3=rows3, rows4=rows4, rows5a=rows5a, rows5b=rows5b,
                precision=precision, paths=paths, hashes=hashes)

def outputs_text(d):
    s1 = longtable("S1",
        r"Country-level CMUR slopes and direction classes, 2010--2022. A is released annual CMUR; B1 uses admissible source anchors. Slopes are in pp/year. Primary is the common two-sided 95\% rule; secondary is the common directional 95\% rule, applied to the same fixed fits. CHANGE denotes resolved-class disagreement; UNRESOLVED is an unmet support criterion, not no change. EU27\_2020 is separate from the Member States.",
        ["Geo", "A slope", "B1 slope", "B1 $n$", r"\shortstack{Primary\\A class}", r"\shortstack{Primary\\B1 class}", r"\shortstack{Primary\\comparison}", r"\shortstack{Secondary\\A class}", r"\shortstack{Secondary\\B1 class}", r"\shortstack{Secondary\\comparison}"],
        "lrrrllllll", d["rows1"], size=r"\scriptsize", padding=3)
    a, b, ratio, wider = d["precision"]
    notes2 = rf"\noindent Member-State summaries from unrounded intervals: median A width {a:.4f} pp/year; median B1 width {b:.4f} pp/year; median paired B1/A ratio {ratio:.3f}. B1 is wider in {wider}/27 Member States; EU27\_2020 is not included in that count."
    s2 = longtable("S2",
        r"Country-level classical model-based OLS 95\% intervals, widths and directional Student-$t$ sign support, 2010--2022. Widths are in pp/year; the ratio is B1/A. Both $p_+$ values describe the fitted model. Neither intervals nor sign support are Bayesian posterior probabilities. Summaries use unrounded values.",
        ["Geo", "A 95\\% interval", "B1 95\\% interval", "A width", "B1 width", "B1/A", "A $p_+$", "B1 $p_+$"],
        "lrrrrrrr", d["rows2"], notes2, size=r"\footnotesize", padding=4)
    s3 = longtable("S3",
        r"Complete tested short-window threshold grid through 2022. U is UNSCALED\_STRESS; B is BROWNIAN\_SD\_1\_SQRT2. A/D/U counts mean agreement/resolved-class disagreement/unresolved; U in the scenario column denotes the unscaled scenario. Shares are percentages of unresolved comparisons, not probabilities of truth. Thresholds are discrete; no interpolation is implied.",
        ["Population", "Metric", "Scenario", r"$\tau$", "A", "D", "U", "Total", r"\shortstack{Unresolved\\share (\%)}"],
        "lllrrrrrr", d["rows3"], size=r"\footnotesize", padding=3)
    s4 = longtable("S4",
        r"Continuous weakest-variant support: conditional, descriptive distributions within each scenario's concordant preferred-class subset. U and B use the same scenario definitions as Table~\ref{tab:S3}; their concordant subsets need not be paired. These statistics are not posterior, joint or calibrated probabilities. Conflicts, ties and missing coverage are excluded from support quantiles and reported separately. For turning, support uses the lower constituent-change support. Quantiles use the concordant denominator $n$.",
        ["Population", "Metric", "Scenario", "Total", r"\shortstack{Concordant\\$n$}", "Conflicts", "Ties", "Missing", "P10", "Q25", "Median", "Q75", "P90"],
        "lllrrrrrrrrrr", d["rows4"], size=r"\footnotesize", padding=3)
    # One numbered longtable with explicit panels; the panel-B values occupy two columns.
    # Build S5 separately so structural multicolumn rows are not escaped.
    headers = [r"Diagnostic / statistic", "Population / value", r"\shortstack{New source-aware\\resolutions}", r"\shortstack{New final\\resolutions}"]
    ordinary = longtable("S5",
        r"Supporting diagnostics through 2022. Panel A uses the unscaled reference scenario at $\tau=0.95$ and counts additional resolved acceleration claims when the deterministic midpoint vote is removed: the literal diagnostic retains minimum three applicable variants; the separately labelled exploratory stochastic-only diagnostic requires at least two. Both require unanimity among applicable votes. Zero additional resolutions show that midpoint geometry alone is insufficient under this design, without identifying a causal share. Panel B pairs the same geography and central observed target, with all required anchors admissible. Multi-span windows are dependent and descriptive; the comparison supplies no scaling exponent, Brownian validation, true two-year calibration or casewise monotonicity.",
        headers, "p{6.1cm}lrr", [], size=r"\footnotesize", padding=3)
    insert = [r"\multicolumn{4}{l}{\textit{Panel A: midpoint attribution}}" + END]
    insert += [line(row) for row in d["rows5a"]]
    insert += [r"\midrule", r"\multicolumn{4}{l}{\textit{Panel B: paired Member-State multi-span diagnostic}}" + END]
    insert += [tex(row[0]) + " & " + tex(row[1]) + r" & -- & --" + END for row in d["rows5b"]]
    s5 = ordinary.replace(r"\endlastfoot" + "\n", r"\endlastfoot" + "\n" + "\n".join(insert) + "\n")
    return (s1, s2, s3, s4, s5)

def main():
    d = prepare()
    texts = outputs_text(d)
    SUP.mkdir(exist_ok=True)
    for path, content in zip(OUTPUTS, texts):
        path.write_text(content, encoding="utf-8", newline="\n")
    require(all(digest(d["paths"][key]) == h for key, h in d["hashes"].items()), "Protected source changed")
    print(json.dumps(dict(source_sha256=d["hashes"], output_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in OUTPUTS},
                         rows=dict(S1=28, S2=28, S3=84, S4=12, S5_midpoint=4, S5_pair_statistics=7),
                         scientific_analysis_rerun=False), indent=2))

if __name__ == "__main__":
    main()
