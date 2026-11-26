# KinaseVarDB — Cancer Kinase Variants & Inhibitors

An end-to-end data pipeline and interactive Streamlit dashboard exploring UniProt protein domains, ClinVar variants, and ChEMBL targeted inhibitors across 8 cancer-relevant kinase genes (`BRAF`, `ALK`, `KIT`, `RET`, `MET`, `PDGFRA`, `ERBB2`, and `EGFR`).

## Project Structure

- `schema.sql` — SQLite schema defining `proteins`, `domains`, `variants`, and `drugs`.
- `data/raw/` — Raw CSV datasets (`proteins.csv`, `domains.csv`, `variants.csv`).
- `data/processed/` — Transformed datasets (`clean_proteins.csv`, `clean_domains.csv`, `clean_variants.csv`).
- `etl/`
  - `01_fetch.py` — Fetches raw protein, domain, and variant data into `data/raw/`.
  - `02_transform.py` — Parses HGVS protein, range, and frameshift positions (`aa_position`), flags kinase-domain overlap (`in_kinase_domain`, `matched_domain`), and writes `data/processed/`.
  - `03_load.py` — Validates columns against SQLite `PRAGMA` metadata, deduplicates records, and loads `kinasevardb.sqlite` idempotently.
  - `04_fetch_drugs.py` — Fetches targeted kinase inhibitors (Phase ≥ 2) and loads `drugs` data.
- `app.py` — Streamlit dashboard with gene selection, kinase-domain filtering, an interactive Plotly 1D Domain & Variant Map, domain/variant/inhibitor tables, and a cross-gene variant summary.
- `kinasevardb.sqlite` — Local SQLite database queried by `app.py`.
- `requirements.txt` — Python dependencies for the pipeline and dashboard.

## Setup & Usage

1. Create and activate a virtual environment, then install dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt