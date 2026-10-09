"""Deterministic, outcome-blind WP1A transformation of frozen JSON-stat slices."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from acquire_wp1a import GEOS, ROOT, SLICES

OUT = ROOT / "data" / "processed"
META = ROOT / "data" / "metadata"


def read_sources():
    manifest = json.loads((ROOT / "logs" / "acquisition_manifest.json").read_text(encoding="utf-8"))
    sources = {}
    for rec in manifest:
        path = ROOT / rec["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == rec["sha256"], path
        if rec["name"] in SLICES:
            sources[rec["name"]] = json.loads(path.read_text(encoding="utf-8"))
    assert set(sources) == set(SLICES)
    return manifest, sources


def cells(data):
    ids, sizes = data["id"], data["size"]
    lookup = []
    for dim in ids:
        index = data["dimension"][dim]["category"]["index"]
        lookup.append({int(v): k for k, v in index.items()})
    values, flags = data.get("value", {}), data.get("status", {})
    result = {}
    n = math.prod(sizes)
    for flat in range(n):
        remainder, coords = flat, []
        for size in reversed(sizes):
            coords.append(remainder % size)
            remainder //= size
        dims = dict(zip(ids, (l[i] for l, i in zip(lookup, reversed(coords)))))
        key = (dims["geo"], int(dims["time"]), dims.get("indic_env", ""))
        assert key not in result, key
        result[key] = (values.get(str(flat)), flags.get(str(flat), ""))
    return result


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def inventory(manifest, sources):
    rows = []
    for ds, data in sources.items():
        for dim in data["id"]:
            cat = data["dimension"][dim]["category"]
            for code, pos in cat["index"].items():
                rows.append(dict(dataset=ds, dimension_id=dim, code=code, position=pos,
                                 label=cat.get("label", {}).get(code, ""), selection=SLICES[ds].get(dim, "all requested geographies / 2010-2024"),
                                 dataset_updated=data.get("updated", "")))
    write_csv(META / "dimension_codelist_inventory.csv", rows,
              ["dataset", "dimension_id", "code", "position", "label", "selection", "dataset_updated"])
    dsd_rows, flag_rows = [], []
    for rec in manifest:
        if not rec["name"].endswith("_dsd"):
            continue
        ds = rec["name"].removesuffix("_dsd")
        root = ET.parse(ROOT / rec["path"]).getroot()
        for item in root.iter():
            tag = item.tag.split("}")[-1]
            if tag in {"Dataflow", "DataStructure", "Codelist", "ConceptScheme"}:
                dsd_rows.append(dict(dataset=ds, object_type=tag, object_id=item.get("id", ""), version=item.get("version", "")))
            if tag == "Codelist" and item.get("id") == "OBS_FLAG":
                for code in item:
                    if code.tag.split("}")[-1] != "Code":
                        continue
                    label = next((x.text for x in code if x.tag.split("}")[-1] == "Name" and x.get("{http://www.w3.org/XML/1998/namespace}lang") == "en"), "")
                    flag_rows.append(dict(dataset=ds, code=code.get("id"), label=label, codelist_version=item.get("version")))
    write_csv(META / "dsd_versions.csv", dsd_rows, ["dataset", "object_type", "object_id", "version"])
    write_csv(META / "observation_flag_definitions.csv", flag_rows, ["dataset", "code", "label", "codelist_version"])
    (META / "unit_conversion.json").write_text(json.dumps({"working_mass_unit": "thousand tonnes",
      "env_wastrt": {"source_unit": "T", "factor": 0.001, "operation": "divide tonnes by 1000"},
      "env_ac_mfa": {"source_unit": "THS_T", "factor": 1},
      "env_ac_sd": {"source_unit": "THS_T", "factor": 1},
      "env_ac_cur": {"source_unit": "PC", "factor": 1, "precision": "one decimal percentage point"}}, indent=2), encoding="utf-8")


def main():
    manifest, src = read_sources()
    inventory(manifest, src)
    lookup = {name: cells(data) for name, data in src.items()}
    panel, gaps, exceptions, anchors, support = [], [], [], [], []
    waste = lookup["env_wastrt"]
    for geo in GEOS:
        for year in range(2010, 2025):
            cm, cmflag = lookup["env_ac_cur"].get((geo, year, ""), (None, ""))
            wr, wf = waste.get((geo, year, ""), (None, ""))
            dm, df = lookup["env_ac_mfa"].get((geo, year, "DMC"), (None, ""))
            imp, imf = lookup["env_ac_sd"].get((geo, year, "IMP_RCV_R"), (None, ""))
            exp, exf = lookup["env_ac_sd"].get((geo, year, "EXP_RCV_R"), (None, ""))
            raw_observed = wr is not None
            admissible = raw_observed and not wf
            if raw_observed and wf:
                provenance = "source_exception"
            elif raw_observed:
                provenance = "waste_anchor_reported"
            elif year >= 2023:
                provenance = "endpoint_unanchored"
            elif year % 2:
                provenance = "regular_intervening_estimated"
            else:
                provenance = "source_exception"
            reasons = []
            if cm is None: reasons.append("released_cmur_missing")
            if dm is None: reasons.append("dmc_missing")
            if imp is None: reasons.append("import_correction_missing")
            if exp is None: reasons.append("export_correction_missing")
            if year % 2 == 0 and not raw_observed: reasons.append("raw_waste_anchor_missing")
            if wf: reasons.append(f"waste_flag_{wf}")
            if cmflag: reasons.append(f"cmur_flag_{cmflag}")
            if df: reasons.append(f"dmc_flag_{df}")
            if imf: reasons.append(f"import_flag_{imf}")
            if exf: reasons.append(f"export_flag_{exf}")
            wm = wr / 1000 if wr is not None else None
            u = wm - imp + exp if None not in (wm, imp, exp) else None
            reconstructed = 100 * u / (dm + u) if u is not None and dm is not None and dm + u > 0 else None
            gap = reconstructed - cm if reconstructed is not None and cm is not None else None
            # Published precision envelope: waste to nearest tonne, DMC to
            # 0.001 thousand tonnes, trade flows to whole thousand tonnes,
            # released CMUR to 0.1 percentage point.
            precision_consistent = None
            precision_excess = None
            if reconstructed is not None and cm is not None:
                u_delta = 0.0005 + 0.5 + 0.5
                low = 100 * (u-u_delta) / (dm+0.0005+u-u_delta)
                high = 100 * (u+u_delta) / (dm-0.0005+u+u_delta)
                precision_excess = max(0.0, low-(cm+0.05), (cm-0.05)-high)
                precision_consistent = precision_excess == 0
            row = dict(geo=geo, year=year, cmur_released_pct=cm, cmur_released_flag=cmflag,
                       waste_rcv_r_mass=wm, waste_source_flag=wf, waste_provenance_class=provenance,
                       waste_anchor_observed=admissible, dmc_mass=dm, dmc_flag=df, impw_mass=imp,
                       expw_mass=exp, trade_source_flag="|".join(x for x in (imf, exf) if x),
                       u_reconstructed_mass=u, cmur_reconstructed_pct=reconstructed,
                       reconstruction_gap_pp=gap,
                       source_aware_status="unresolved_linkage" if admissible and None not in (cm, dm, imp, exp) else "unresolved_source",
                       exclusion_reason=";".join(reasons))
            panel.append(row)
            anchors.append(dict(geo=geo, year=year, raw_waste_value_tonnes=wr, raw_waste_flag=wf,
                                raw_observation_present=raw_observed, admissible_anchor=admissible,
                                provenance_class=provenance))
            if gap is not None:
                gaps.append(dict(geo=geo, year=year, released_pct=cm, reconstructed_pct=reconstructed,
                                 gap_pp=gap, precision_consistent=precision_consistent,
                                 precision_excess_pp=precision_excess, waste_flag=wf, cmur_flag=cmflag, dmc_flag=df,
                                 import_flag=imf, export_flag=exf))
            if reasons or provenance != "waste_anchor_reported":
                exceptions.append(dict(geo=geo, year=year, provenance_class=provenance,
                                       raw_anchor_present=raw_observed, reasons=";".join(reasons)))
    for geo in GEOS:
        for year in range(2012, 2021, 2):
            trio = [waste.get((geo, t, ""), (None, "")) for t in (year-2, year, year+2)]
            reasons = []
            for t, (value, flag) in zip((year-2, year, year+2), trio):
                if value is None: reasons.append(f"{t}_missing")
                if flag: reasons.append(f"{t}_flag_{flag}")
            if not reasons and (trio[0][0] + trio[2][0]) / 2 <= 0:
                reasons.append("nonpositive_neighbor_mean")
            support.append(dict(geo=geo, interior_anchor_year=year, eligible=not reasons,
                                ineligible_reason=";".join(reasons)))
    count = Counter(r["geo"] for r in support if r["eligible"])
    counts = [dict(geo=g, eligible_residual_count=count[g], possible_interior_anchors=5,
                   ineligible_reasons=";".join(f'{r["interior_anchor_year"]}:{r["ineligible_reason"]}' for r in support if r["geo"] == g and not r["eligible"])) for g in GEOS]
    counts.append(dict(geo="POOLED_EU27_MS", eligible_residual_count=sum(count[g] for g in GEOS if g != "EU27_2020"),
                       possible_interior_anchors=135, ineligible_reasons=""))
    write_csv(OUT / "cmur_source_provenance_panel.csv", panel, list(panel[0]))
    write_csv(OUT / "verified_anchor_map.csv", anchors, list(anchors[0]))
    write_csv(OUT / "reconstruction_gap.csv", gaps, list(gaps[0]))
    write_csv(OUT / "exclusion_flag_register.csv", exceptions, list(exceptions[0]))
    write_csv(OUT / "pseudogap_eligibility.csv", support, list(support[0]))
    write_csv(OUT / "pseudogap_support.csv", counts, list(counts[0]))
    summary = {"panel_rows": len(panel), "raw_waste_cells": sum(r["raw_observation_present"] for r in anchors),
               "admissible_anchors": sum(r["admissible_anchor"] for r in anchors),
               "reconstruction_cells": len(gaps), "pooled_pseudogap_count": counts[-1]["eligible_residual_count"],
               "absolute_gap_pp_max": max(abs(r["gap_pp"]) for r in gaps),
               "absolute_gap_pp_median": sorted(abs(r["gap_pp"]) for r in gaps)[len(gaps)//2],
               "gap_over_0_05_pp": sum(abs(r["gap_pp"]) > 0.05 for r in gaps),
               "precision_inconsistent_cells": sum(not r["precision_consistent"] for r in gaps),
               "flagged_waste_cells": [(r["geo"], r["year"], r["raw_waste_flag"]) for r in anchors if r["raw_waste_flag"]]}
    (ROOT / "logs" / "transformation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
