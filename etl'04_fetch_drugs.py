import time
import requests

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
def fetch_drugs_for_gene(gene: str, target_chembl_id: str) -> list[dict]:
    resp = requests.get(
        f"{CHEMBL_BASE}/mechanism.json",
        params={"target_chembl_id": target_chembl_id, "limit": 100},
        timeout=30,
    )
    resp.raise_for_status()
    mechanisms = resp.json().get("mechanisms", [])

    rows = []
    seen_mols = set()

    for mech in mechanisms:
        mol_id = mech.get("molecule_chembl_id")
        if not mol_id or mol_id in seen_mols:
            continue

        m_resp = requests.get(f"{CHEMBL_BASE}/molecule/{mol_id}.json", timeout=30)
        if m_resp.status_code != 200:
            continue
        mol = m_resp.json()

        raw_phase = mol.get("max_phase")
        try:
            max_phase = float(raw_phase) if raw_phase is not None else 0.0
        except (TypeError, ValueError):
            max_phase = 0.0

        # Keep approved (phase 4) or clinical-stage (phase >= 2) drugs
        if max_phase >= 2:
            seen_mols.add(mol_id)
            pref_name = mol.get("pref_name") or mol_id
            rows.append({
                "gene_symbol": gene,
                "chembl_molecule_id": mol_id,
                "drug_name": pref_name.title(),
                "max_phase": int(max_phase),
                "mechanism_of_action": mech.get("mechanism_of_action", ""),
            })
        time.sleep(0.1)

    return rows
def save_drugs_to_sqlite(rows: list[dict], db_path: str = "kinasevardb.sqlite") -> None:
    import sqlite3
    import pandas as pd

    df = pd.DataFrame(rows)
    conn = sqlite3.connect(db_path)
    try:
        if not df.empty:
            df.to_sql("drugs", conn, if_exists="replace", index=False)
        else:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS drugs (
                    gene_symbol TEXT,
                    chembl_molecule_id TEXT,
                    drug_name TEXT,
                    max_phase INTEGER,
                    mechanism_of_action TEXT
                )
                """
            )
        conn.commit()
    finally:
        conn.close()


def main() -> None:
    fallback_targets = {
        "EGFR": "CHEMBL203",
        "BRAF": "CHEMBL5145",
        "ALK": "CHEMBL4247",
    }
    genes = globals().get("GENES", ["EGFR", "BRAF", "ALK"])
    all_rows = []

    for gene in genes:
        print(f"Fetching ChEMBL drugs for {gene}...")
        if "get_target_chembl_id" in globals():
            target_id = get_target_chembl_id(gene)
        else:
            target_id = fallback_targets.get(gene)

        if not target_id:
            print(f"  No ChEMBL target ID found for {gene}")
            continue

        gene_drugs = fetch_drugs_for_gene(gene, target_id)
        print(f"  Found {len(gene_drugs)} drugs (phase >= 2) for {gene}")
        all_rows.extend(gene_drugs)

    save_drugs_to_sqlite(all_rows)
    print(f"Saved {len(all_rows)} total drug rows to 'drugs' table in kinasevardb.sqlite.")


if __name__ == "__main__":
    main()
