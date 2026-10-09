# CMUR Temporal-Resolution Reproducibility Package

This repository contains archived statistical inputs, provenance data, analysis code, verification scripts, tables and figures for a temporal-resolution analysis of the EU Circular Material Use Rate (CMUR). It covers the 27 EU Member States over **2010–2022**, with **EU27_2020** analysed separately. The source-acquisition and provenance panel span **2010–2024**.

**Creator:** Krisztián Attila Bakon  
**Affiliation:** University of Pannonia – Nagykanizsa University Center for Circular Economy  
**ORCID:** https://orcid.org/0000-0003-1882-094X  
**Contact:** bakon.krisztian@pen.uni-pannon.hu

## Repository structure

| Path | Content |
|---|---|
| `data/raw/` | Four filtered Eurostat JSON-stat snapshots, acquired 26 September 2026. |
| `data/metadata/` | Four Eurostat SDMX structural snapshots plus derived codelist, flag and unit extracts. |
| `data/processed/` | Provenance panel, anchor map, analysis tables, residuals and classifications. |
| `scripts/` | Acquisition, transformation, analysis and independent verification code. |
| `manuscript/` | **Technical asset path retained by the figure/table builder**: generated figures, LaTeX tables and their source code. |
| `logs/` | Filtered source-acquisition manifest and analysis verification state. |
| `SOURCE_CATALOG.csv` | Official source URLs, query filters, UTC retrieval times, archived hashes and rights treatment for all 13 source records. |
| `MANIFEST.csv`, `CHECKSUMS.sha256` | File roles, origins, rights statuses and SHA-256 integrity checks. |

## Sources and rights

The archived Eurostat source slices are unchanged official data/metadata. The processed tables, figures, code and documentation are author-generated. Eurostat permits reuse of statistical data and metadata with attribution, subject to its published exceptions: https://ec.europa.eu/eurostat/help/copyright-notice . Source: Eurostat, datasets `env_ac_cur`, `env_wastrt`, `env_ac_mfa` and `env_ac_sd`, filtered as recorded in `logs/acquisition_manifest.json`; accessed **26 September 2026**. The derived material adapts/analyzes those inputs; Eurostat is not responsible for the interpretations.

Eight approved Eurostat statistical/structural snapshots are included. Four ancillary ESMS HTML pages are linked rather than redistributed. The 2018 Eurostat CMUR calculation-method PDF is linked but excluded because the full publication credits a third-party photograph. `SOURCE_CATALOG.csv` gives the original URLs and archived study hashes. Author-written code is MIT licensed; author-generated derived data, documentation, tables and figures are CC BY 4.0 licensed. Official Eurostat snapshots retain Eurostat's separate reuse conditions. See `LICENSES.md`, `LICENSE-CODE-MIT.txt`, `LICENSE-DATA-DOCS-CC-BY-4.0.txt` and `THIRD_PARTY_NOTICES.md`. No European Environment Agency raw file is included.

## Runtime

The byte-identical clean-room analysis was tested on 1 October 2026 with Python **3.13.15**, NumPy **2.2.6**, SciPy **1.16.0** and ReportLab **4.4.9**. Install the tested dependencies with `python -m pip install -r requirements.txt`. Exact interpreter/library versions from the initial analysis were not archived. A Python 3.12 / NumPy 1.26 / SciPy 1.13 test retained headline classes but changed some last-digit slope serialization, so use the pinned environment for exact file hashes. Poppler `pdftoppm` **26.07.0** is optional for PNG previews; vector PDFs and LaTeX tables require ReportLab but not Poppler. No R runtime is used.

## Reproduce the archived analysis

Run these commands from the repository root. They write only inside a local copy of this tree. The temporal simulation uses random seed **20260926** and **20,000** trajectories per applicable stochastic variant.

```text
python -m pip install -r requirements.txt
python scripts/verify_checksums.py
python scripts/build_wp1a.py
python scripts/verify_wp1a.py
python scripts/wp1b_r.py --gate-only
python scripts/wp1b_r.py --analyze
python scripts/verify_wp1b_r.py
python scripts/check_outputs.py
python manuscript/scripts/build_assets.py
```

Run `verify_checksums.py` **before** regenerating outputs. The source and temporal verifiers check structure, arithmetic and deterministic reruns; `check_outputs.py` checks the numerical summaries. The asset builder recreates the included figures and tables from the processed CSVs. It may overwrite assets in the local copy. `MANIFEST.csv` records the original hashes for comparison.

Expected checks include **420** geography-year provenance rows for 2010–2024, **193** admissible raw waste anchors, **130** pooled Member-State interior residuals, medium-term Member-State classes **15 agreement / 0 resolved-class disagreement / 12 unresolved**, annual-change classes **10 / 0 / 314**, turning **2 / 0 / 295**, and acceleration **1 / 0 / 296**. EU27_2020 annual-series and anchor-only slopes round to **0.0577** and **0.0536** percentage points/year. These are verification targets, not values substituted into the analysis.

## Reacquiring official sources

The included snapshots make the archived analysis reproducible without a live API. `SOURCE_CATALOG.csv` and `logs/acquisition_manifest.json` preserve exact URLs, filters, UTC retrieval times, dataset updates where available and source SHA-256 hashes. To fetch current official versions, use the listed URLs or run `python scripts/acquire_wp1a.py` **in a separate copy**. Eurostat's public API serves current data, so a fresh response may differ from the archived vintage. Do not overwrite the included snapshots when reproducing the reported run.

The source snapshots support provenance reconstruction, and the author-generated panel/anchor map support the temporal analysis. This package does not claim to reconstruct the released indicator from same-vintage components or to provide an exact annual CMUR error distribution.
