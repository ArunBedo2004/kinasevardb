# KinaseVarDB

An end-to-end bioinformatics data pipeline and interactive Streamlit dashboard exploring cancer kinase protein domains, ClinVar variants, and ChEMBL targeted inhibitors.

- **Live App:** https://kinasevardb.streamlit.app
- **GitHub Repository:** https://github.com/ArunBedo2004/kinasevardb

## Overview

**KinaseVarDB** integrates data from **UniProt**, **ClinVar**, and **ChEMBL** into a SQLite database (`kinasevardb.sqlite`) and provides an interactive web interface (`app.py`) to inspect protein domain architecture, amino acid variant positions, kinase-domain overlap, and targeted drugs across cancer kinases.

### Key Features
- **Multi-Source ETL Pipeline:** Numbered scripts under `etl/` (`01_fetch.py`, `02_transform.py`, `03_load.py`, and `04_fetch_drugs.py`) that fetch, transform, load, and enrich kinase domain, variant, and drug data.
- **Relational SQLite Backend:** Defined in `schema.sql` and stored in `kinasevardb.sqlite` for fast querying.
- **Interactive Dashboard (`app.py`):** Built with Streamlit and Plotly to switch kinase targets, filter variants within the kinase domain, visualize 1D domain-variant maps, and inspect domain, variant, and targeted inhibitor tables.

## Project Structure

- `app.py` — Streamlit and Plotly web application
- `kinasevardb.sqlite` — SQLite database used by the app
- `schema.sql` — Database table definitions
- `requirements.txt` — Python dependencies (`requests`, `pandas`, `streamlit`, `plotly`)
- `README.md` — Project documentation
- `.gitignore` — Git ignore rules
- `etl/`
  - `01_fetch.py` — Fetches source data
  - `02_transform.py` — Cleans and transforms records
  - `03_load.py` — Loads data into SQLite
  - `04_fetch_drugs.py` — Fetches targeted inhibitor data
- `data/`
  - `raw/` — Raw CSV files (`domains.csv`, `proteins.csv`, `variants.csv`)
  - `processed/` — Processed datasets

## Running Locally

1. Clone the repository:
   ```bash
   git clone [https://github.com/ArunBedo2004/kinasevardb.git](https://github.com/ArunBedo2004/kinasevardb.git)
   cd kinasevardb
