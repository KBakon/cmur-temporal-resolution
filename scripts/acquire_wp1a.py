"""Freeze filtered Eurostat WP1A responses, structural metadata and checksums."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
META = ROOT / "data" / "metadata"
LOG = ROOT / "logs"
GEOS = "AT BE BG HR CY CZ DK EE FI FR DE EL HU IE IT LV LT LU MT NL PL PT RO SK SI ES SE EU27_2020".split()
BASE = "https://ec.europa.eu/eurostat/api/dissemination"
WINDOW = {"lang": "EN", "sinceTimePeriod": "2010", "untilTimePeriod": "2024", "geo": GEOS}
SLICES = {
    "env_ac_cur": {"freq": "A", "unit": "PC"},
    "env_wastrt": {"freq": "A", "unit": "T", "hazard": "HAZ_NHAZ", "wst_oper": "RCV_R", "waste": "TOTAL"},
    "env_ac_mfa": {"freq": "A", "indic_env": "DMC", "material": "TOTAL", "unit": "THS_T"},
    "env_ac_sd": {"freq": "A", "indic_env": ["IMP_RCV_R", "EXP_RCV_R"], "material": "TOTAL", "unit": "THS_T"},
}


def freeze(name: str, url: str, suffix: str, directory: Path) -> dict:
    req = Request(url, headers={"User-Agent": "CMUR-TEMPORAL-AUDIT-WP1A/1.0"})
    with urlopen(req, timeout=180) as response:
        body = response.read()
        headers = {k: v for k, v in response.headers.items() if k.lower() in {"last-modified", "etag", "content-type"}}
    sha = hashlib.sha256(body).hexdigest()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}_{sha[:16]}.{suffix}"
    if path.exists():
        if path.read_bytes() != body:
            raise RuntimeError(f"Hash collision: {path}")
    else:
        path.write_bytes(body)
    record = {"name": name, "url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
              "sha256": sha, "bytes": len(body), "path": str(path.relative_to(ROOT)).replace('\\', '/'),
              "http_headers": headers}
    if suffix == "json":
        data = json.loads(body)
        record.update(dataset_updated=data.get("updated"), dimension_ids=data.get("id"), dimension_sizes=data.get("size"))
    return record


def main() -> None:
    records = []
    for ds, selected in SLICES.items():
        query = dict(WINDOW, **selected)
        url = f"{BASE}/statistics/1.0/data/{ds}?{urlencode(query, doseq=True)}"
        record = freeze(ds, url, "json", RAW)
        records.append(record)
        print(ds, record["bytes"], record.get("dataset_updated"), flush=True)
        records.append(freeze(f"{ds}_dsd", f"{BASE}/sdmx/2.1/dataflow/ESTAT/{ds}?references=descendants", "xml", META))
        metadata_name = {"env_wastrt": "env_wasgt_esms", "env_ac_mfa": "env_ac_mfa_sims"}.get(ds, f"{ds}_esms")
        records.append(freeze(f"{ds}_esms", f"https://ec.europa.eu/eurostat/cache/metadata/en/{metadata_name}.htm", "html", META))
    records.append(freeze("cmur_calculation_method_2018", "https://ec.europa.eu/eurostat/documents/3859598/9407565/KS-FT-18-009-EN-N.pdf", "pdf", META))
    LOG.mkdir(exist_ok=True)
    path = LOG / "acquisition_manifest.json"
    if path.exists():
        old = path.read_bytes()
        if json.loads(old) != records:
            archive = LOG / f"acquisition_manifest_prior_{hashlib.sha256(old).hexdigest()[:16]}.json"
            if not archive.exists():
                archive.write_bytes(old)
    path.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
