"""Read-only checks of the deposited CMUR analysis products."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from wp1b_r import N_DRAWS, SEED


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def expect(label: str, got: object, wanted: object) -> None:
    if got != wanted:
        raise AssertionError(f"{label}: got {got!r}; expected {wanted!r}")
    print(f"PASS {label}: {got}")


def main() -> None:
    panel = rows("cmur_source_provenance_panel.csv")
    anchor = rows("verified_anchor_map.csv")
    residual = rows("b2r_residual_pool.csv")
    comparison = rows("temporal_stability_comparison.csv")
    a = {r["geo"]: r for r in rows("branchA_medium_term.csv")}
    b = {r["geo"]: r for r in rows("branchB1_anchor_medium_term.csv")}
    expect("panel rows 2010–2024", len(panel), 420)
    expect("admissible raw waste anchors", sum(r["admissible_anchor"] == "True" for r in anchor), 193)
    expect("pooled Member-State residuals", sum(r["pooled_member_state"] == "True" for r in residual), 130)
    expected = {
        "medium_term_slope": (2010, 2022, {"AGREEMENT": 15, "CHANGE": 0, "UNRESOLVED": 12}),
        "annual_change": (2011, 2022, {"AGREEMENT": 10, "CHANGE": 0, "UNRESOLVED": 314}),
        "turning": (2012, 2022, {"AGREEMENT": 2, "CHANGE": 0, "UNRESOLVED": 295}),
        "acceleration": (2012, 2022, {"AGREEMENT": 1, "CHANGE": 0, "UNRESOLVED": 296}),
    }
    for metric, (start, end, target) in expected.items():
        observed = Counter(
            r["stability_class"] for r in comparison
            if r["geo"] != "EU27_2020" and r["metric"] == metric
            and (metric == "medium_term_slope" or start <= int(r["year"]) <= end)
        )
        expect(metric + " Member-State classes", {k: observed[k] for k in target}, target)
    expect("EU27_2020 annual slope, 4 decimals", round(float(a["EU27_2020"]["slope_pp_per_year"]), 4), 0.0577)
    expect("EU27_2020 anchor slope, 4 decimals", round(float(b["EU27_2020"]["slope_pp_per_year"]), 4), 0.0536)
    expect("simulation seed", SEED, 20260926)
    expect("stochastic trajectories per variant", N_DRAWS, 20000)


if __name__ == "__main__":
    main()
