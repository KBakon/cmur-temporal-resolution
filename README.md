# CMUR Temporal-Resolution Reproducibility Package

Version: `v1.1.0`

Science version: `CMUR-TEMPORAL-AUDIT-P2-2026-11-11`

This repository contains the frozen statistical inputs, provenance materials, processed analysis tables, analysis and verification scripts, diagnostic evidence, and manuscript figure/table assets supporting the CMUR temporal-resolution study.

**Creator:** Krisztián Attila Bakon
**Affiliation:** University of Pannonia – Nagykanizsa University Center for Circular Economy
**ORCID:** https://orcid.org/0000-0003-1882-094X
**Contact:** bakon.krisztian@pen.uni-pannon.hu

## Repository structure

| Path | Content |
|---|---|
| `data/raw/` | Archived filtered Eurostat statistical snapshots |
| `data/metadata/` | Eurostat structural metadata and derived metadata tables |
| `data/processed/` | Frozen provenance, anchor, residual and temporal-analysis tables |
| `analysis/wp_p1a_temporal_support_diagnostics/` | Frozen temporal-support diagnostic evidence and verification |
| `analysis/wp_p1c_calibration_robustness/` | Frozen calibration and robustness diagnostic evidence and verification |
| `scripts/` | Acquisition, transformation, analysis and verification code |
| `manuscript/` | Manuscript source sections, bibliography, figures, tables and supplementary assets |
| `logs/` | Frozen source-acquisition manifest |
| `provenance/SOURCE_CATALOG.csv` | Official source URLs, retrieval information, hashes and rights treatment |
| `MANIFEST_SHA256.txt` | File-level SHA-256 integrity manifest |

## Reproduction environment

The archived `requirements.txt` records the tested NumPy, SciPy and ReportLab versions.

Install dependencies in a disposable environment:

    python -m pip install -r requirements.txt

The frozen repository includes the production scripts and two saved diagnostic workflows.

P1A diagnostic workflow:

    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py baseline
    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py diagnostics
    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py verify

P1C diagnostic workflow:

    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py baseline
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py run
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py report

Production and manuscript-asset entry points include:

    scripts/build_wp1a.py
    scripts/wp1b_r.py
    scripts/verify_wp1a.py
    scripts/verify_wp1b_r.py
    manuscript/scripts/build_assets.py
    manuscript/scripts/build_supplement_assets.py
    manuscript/scripts/verify_integration.py
    manuscript/scripts/verify_supplement_integration.py

Run reproduction or diagnostic commands only in a disposable checkout. Do not overwrite the archived frozen inputs when testing current upstream data.

## Data provenance and rights

The included raw statistical and structural snapshots are official Eurostat files for `env_ac_cur`, `env_wastrt`, `env_ac_mfa`, and `env_ac_sd`.

Exact source URLs, retrieval times, filters and hashes are recorded in `provenance/SOURCE_CATALOG.csv` and `logs/acquisition_manifest.json`.

Eurostat material remains subject to Eurostat reuse conditions and attribution. Author-written code is released under the MIT License. Author-generated derived data, documentation, figures and tables are released under CC BY 4.0.

The repository therefore has mixed rights; see `LICENSES.md`, `LICENSE-CODE-MIT.txt`, `LICENSE-DATA-DOCS-CC-BY-4.0.txt`, and `THIRD_PARTY_NOTICES.md`.

No EEA raw data file is included.

## Scope and boundary

This release preserves the frozen inputs, scripts, diagnostics, verification evidence, and table/figure assets used by the submission-stage CMUR analysis.

It does not reconstruct same-vintage CMUR from components, identify the true intervening-year error distribution, or add empirical validation beyond the documented frozen analyses.

## Release history and citation

The initial `v1.0.0` release is archived on Zenodo:

https://doi.org/10.5281/zenodo.23263064

The present `v1.1.0` release extends the initial public package with frozen P1A and P1C diagnostics, expanded verification outputs, manuscript source materials, and supplementary Tables S1–S5.

The version-specific Zenodo DOI for `v1.1.0` will be added after the GitHub release has been archived by Zenodo.

The `v1.0.0` DOI identifies only the initial release and should not be used to cite the expanded `v1.1.0` contents.
