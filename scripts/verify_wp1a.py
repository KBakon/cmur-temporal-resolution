"""Fail-loud WP1A structure, provenance, arithmetic and rerun verification."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from acquire_wp1a import GEOS, ROOT, SLICES
from build_wp1a import cells, main as rebuild, read_sources

OUT = ROOT / "data" / "processed"
PRODUCTS = ["cmur_source_provenance_panel.csv", "verified_anchor_map.csv", "reconstruction_gap.csv",
            "exclusion_flag_register.csv", "pseudogap_eligibility.csv", "pseudogap_support.csv"]


def table(name):
    with (OUT / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(name, condition, failures, results):
    results[name] = bool(condition)
    if not condition:
        failures.append(name)


def main():
    failures, results = [], {}
    manifest, src = read_sources()
    panel = table(PRODUCTS[0]); anchors = table(PRODUCTS[1]); gaps = table(PRODUCTS[2])
    eligible = table(PRODUCTS[4]); support = table(PRODUCTS[5])
    expected = {(g, y) for g in GEOS for y in range(2010, 2025)}
    keys = [(r["geo"], int(r["year"])) for r in panel]
    check("panel_unique_keys", len(keys) == len(set(keys)), failures, results)
    check("geography_and_time_domain", set(keys) == expected and len(panel) == 420, failures, results)
    check("anchor_map_same_keys", {(r["geo"], int(r["year"])) for r in anchors} == expected, failures, results)
    for ds, selected in SLICES.items():
        data = src[ds]
        check(f"{ds}_dimension_order", data["id"] == ({"env_ac_cur":["freq","unit","geo","time"],
              "env_wastrt":["freq","unit","hazard","wst_oper","waste","geo","time"],
              "env_ac_mfa":["freq","indic_env","material","unit","geo","time"],
              "env_ac_sd":["freq","indic_env","material","unit","geo","time"]}[ds]), failures, results)
        for dim, code in selected.items():
            found = set(data["dimension"][dim]["category"]["index"])
            check(f"{ds}_{dim}_selection", found == set(code if isinstance(code, list) else [code]), failures, results)
        check(f"{ds}_geographies", set(data["dimension"]["geo"]["category"]["index"]) == set(GEOS), failures, results)
        times = set(map(int, data["dimension"]["time"]["category"]["index"]))
        check(f"{ds}_time_domain", times == (set(range(2010, 2023, 2)) if ds == "env_wastrt" else set(range(2010, 2025))), failures, results)
        check(f"{ds}_nonnegative_raw_values", all(v >= 0 for v in data.get("value", {}).values()), failures, results)
    look = {ds: cells(data) for ds, data in src.items()}
    anchor_lookup = {(r["geo"], int(r["year"])): r for r in anchors}
    panel_lookup = {(r["geo"], int(r["year"])): r for r in panel}
    flags_ok = anchors_ok = formula_ok = missing_ok = provenance_ok = True
    for geo, year in expected:
        p = panel_lookup[(geo, year)]; a = anchor_lookup[(geo, year)]
        waste, wf = look["env_wastrt"].get((geo, year, ""), (None, ""))
        cm, cf = look["env_ac_cur"].get((geo, year, ""), (None, ""))
        dm, df = look["env_ac_mfa"].get((geo, year, "DMC"), (None, ""))
        im, imf = look["env_ac_sd"].get((geo, year, "IMP_RCV_R"), (None, ""))
        ex, exf = look["env_ac_sd"].get((geo, year, "EXP_RCV_R"), (None, ""))
        flags_ok &= p["waste_source_flag"] == wf and p["cmur_released_flag"] == cf and p["dmc_flag"] == df and p["trade_source_flag"] == "|".join(x for x in (imf, exf) if x)
        anchors_ok &= (a["raw_observation_present"] == str(waste is not None)) and (a["admissible_anchor"] == str(waste is not None and not wf))
        anchors_ok &= (p["waste_anchor_observed"] == a["admissible_anchor"])
        provenance_ok &= p["waste_provenance_class"] == ("source_exception" if waste is not None and wf else "waste_anchor_reported" if waste is not None else "endpoint_unanchored" if year >= 2023 else "regular_intervening_estimated" if year % 2 else "source_exception")
        if waste is None:
            missing_ok &= p["waste_rcv_r_mass"] == "" and p["u_reconstructed_mass"] == "" and p["cmur_reconstructed_pct"] == ""
        else:
            w = waste / 1000
            u = w-im+ex
            v = 100*u/(dm+u)
            formula_ok &= math.isclose(float(p["waste_rcv_r_mass"]), w, abs_tol=1e-9)
            formula_ok &= math.isclose(float(p["u_reconstructed_mass"]), u, abs_tol=1e-9)
            formula_ok &= math.isclose(float(p["cmur_reconstructed_pct"]), v, abs_tol=1e-12)
            formula_ok &= math.isclose(float(p["reconstruction_gap_pp"]), v-cm, abs_tol=1e-12)
    check("source_flags_preserved", flags_ok, failures, results)
    check("anchors_require_raw_unflagged_observation", anchors_ok, failures, results)
    check("provenance_consistency", provenance_ok, failures, results)
    check("formula_reproducible", formula_ok, failures, results)
    check("missing_waste_not_coerced_to_zero", missing_ok, failures, results)
    check("anchor_coverage", sum(a["raw_observation_present"] == "True" for a in anchors) == 196 and sum(a["admissible_anchor"] == "True" for a in anchors) == 193, failures, results)
    check("gap_rows_are_raw_anchor_cells", len(gaps) == 196 and {(r["geo"],int(r["year"])) for r in gaps} == {(r["geo"],int(r["year"])) for r in anchors if r["raw_observation_present"] == "True"}, failures, results)
    precision_ok = True
    for row in gaps:
        p = panel_lookup[(row["geo"], int(row["year"]))]
        u, dm, cm = float(p["u_reconstructed_mass"]), float(p["dmc_mass"]), float(p["cmur_released_pct"])
        low = 100*(u-1.0005)/(dm+0.0005+u-1.0005)
        high = 100*(u+1.0005)/(dm-0.0005+u+1.0005)
        excess = max(0.0, low-(cm+0.05), (cm-0.05)-high)
        precision_ok &= math.isclose(float(row["precision_excess_pp"]), excess, abs_tol=1e-12)
        precision_ok &= row["precision_consistent"] == str(excess == 0)
    check("published_precision_envelope_reproducible", precision_ok, failures, results)
    counts = Counter(r["geo"] for r in eligible if r["eligible"] == "True")
    support_ok = len(eligible) == 140 and all(int(r["eligible_residual_count"]) == counts[r["geo"]] for r in support if r["geo"] != "POOLED_EU27_MS")
    support_ok &= int(support[-1]["eligible_residual_count"]) == sum(counts[g] for g in GEOS if g != "EU27_2020")
    for row in eligible:
        geo, year = row["geo"], int(row["interior_anchor_year"])
        trio = [look["env_wastrt"].get((geo, y, ""), (None, "")) for y in (year-2, year, year+2)]
        should_qualify = all(v is not None and not flag for v, flag in trio)
        should_qualify &= (trio[0][0] + trio[2][0]) / 2 > 0 if all(v is not None for v, _ in trio) else False
        support_ok &= (row["eligible"] == str(should_qualify))
    check("pseudogap_support_counts", support_ok, failures, results)
    prior = {name: digest(OUT/name) for name in PRODUCTS}
    rebuild()
    after = {name: digest(OUT/name) for name in PRODUCTS}
    check("deterministic_rerun", prior == after, failures, results)
    report = {"status": "PASS" if not failures else "FAIL", "checks": results, "failures": failures,
              "source_sha256": {r["name"]: r["sha256"] for r in manifest},
              "processed_sha256": after, "gap_summary": json.loads((ROOT/"logs"/"transformation_summary.json").read_text(encoding="utf-8"))}
    (ROOT / "logs" / "verification_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(report["status"], len(results), "checks", failures)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
