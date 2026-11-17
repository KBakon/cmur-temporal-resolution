# CMUR Temporal-Resolution Reproducibility Package

Version: `v1.1.0`

Science version: `CMUR-TEMPORAL-AUDIT-P2-2026-11-11`

This repository contains the frozen statistical inputs, provenance materials, processed analysis tables, analysis and verification scripts, diagnostic evidence, supporting the CMUR temporal-resolution study.

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
| `logs/` | Frozen source-acquisition manifest |
| `provenance/SOURCE_CATALOG.csv` | Official source URLs, retrieval information, hashes and rights treatment |
| `MANIFEST_SHA256.txt` | File-level SHA-256 integrity manifest |

## Reproduction environment

The archived `requirements.txt` records the tested NumPy, SciPy and ReportLab versions.

Install dependencies in a disposable environment:

    python -m pip install -r requirements.txt

Run the following from the repository root in a disposable copy of this release (no additional project files are required):

    python -B scripts/verify_wp1a.py
    python -B scripts/verify_wp1b_r.py
    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py baseline
    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py diagnostics
    python -B analysis/wp_p1a_temporal_support_diagnostics/run_diagnostics.py verify
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py baseline
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py run
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py aggregate
    python -B analysis/wp_p1c_calibration_robustness/run_diagnostics.py verify

The production verifiers regenerate processed outputs, so use a disposable copy. Diagnostic baselines run those verifiers again in isolated nested copies and require byte-identical frozen production CSVs. The diagnostic stages reproduce the frozen scientific tables and verify headline gates, paired residual draws, threshold classes, continuous supports and multi-span results. Integrity inventories cover distributed release inputs only, excluding each runner's outputs, mutable verification metadata and the manifest itself to avoid circular hashes. Run P1A before P1C. Source verification follows the public inclusion catalog: the four statistical JSON and four structural XML snapshots are required and hash-checked; the four explicitly excluded ancillary HTML pages are referenced only. The optional P1C `report` stage also redraws diagnostic figures; it is not required for the scientific CSV reproduction above.

## Data provenance and rights

The included raw statistical and structural snapshots are official Eurostat files for `env_ac_cur`, `env_wastrt`, `env_ac_mfa`, and `env_ac_sd`.

Exact source URLs, retrieval times, filters and hashes are recorded in `provenance/SOURCE_CATALOG.csv` and `logs/acquisition_manifest.json`.

Eurostat material remains subject to Eurostat reuse conditions and attribution. Author-written code is released under the MIT License. Author-generated derived data, documentation, figures and tables are released under CC BY 4.0.

The repository therefore has mixed rights; see `LICENSES.md`, `LICENSE-CODE-MIT.txt`, `LICENSE-DATA-DOCS-CC-BY-4.0.txt`, and `THIRD_PARTY_NOTICES.md`.

No EEA raw data file is included.

## Scope and boundary

This release preserves the frozen inputs, scripts, diagnostics, verification evidence, used for the CMUR analysis.

It does not reconstruct same-vintage CMUR from components, identify the true intervening-year error distribution, or add empirical validation beyond the documented frozen analyses.

## Release history and citation

Current version `v1.1.0`: https://doi.org/10.5281/zenodo.24200978

Historical version `v1.0.0`: https://doi.org/10.5281/zenodo.23263064

These version-specific DOIs refer to different releases; cite the current DOI for v1.1.0.
