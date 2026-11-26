from __future__ import annotations

from pathlib import Path
import time
import requests
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

GENES = ["EGFR", "BRAF", "ALK", "ERBB2", "MET", "RET", "KIT", "PDGFRA"]

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "KinaseVarDB/1.0 (bioinformatics-etl)"})


def get_uniprot_data(genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    protein_rows: list[dict] = []
    domain_rows: list[dict] = []

    for gene in genes:
        params = {
            "query": f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)",
            "format": "json",
            "size": 1,
        }
        resp = SESSION.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
        resp.raise_for_status()
        results = resp.json().get("results", [])

        if not results:
            print(f"Warning: no UniProt record found for {gene}")
            continue

        entry = results[0]
        uniprot_id = entry.get("primaryAccession", "")
        seq_len = entry.get("sequence", {}).get("length")
        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", gene)
        )

        protein_rows.append(
            {
                "gene_symbol": gene,
                "uniprot_id": uniprot_id,
                "protein_name": protein_name,
                "sequence_length": seq_len,
            }
        )

        for feat in entry.get("features", []):
            if feat.get("type") in {"Domain", "Region"}:
                desc = feat.get("description", "")
                loc = feat.get("location", {})
                start_val = loc.get("start", {}).get("value")
                end_val = loc.get("end", {}).get("value")
                if start_val is not None and end_val is not None:
                    domain_rows.append(
                        {
                            "gene_symbol": gene,
                            "uniprot_id": uniprot_id,
                            "domain_name": desc or feat.get("type", "Domain"),
                            "start_pos": int(start_val),
                            "end_pos": int(end_val),
                        }
                    )

        time.sleep(0.2)

    return pd.DataFrame(protein_rows), pd.DataFrame(domain_rows)


def _extract_clinical_significance(doc: dict) -> str:
    germline = doc.get("germline_classification", {}) or {}
    if germline.get("description"):
        return str(germline["description"])
    clin_sig = doc.get("clinical_significance", {}) or {}
    if clin_sig.get("description"):
        return str(clin_sig["description"])
    return ""


def get_clinvar_variants(genes: list[str], max_per_gene: int = 150) -> pd.DataFrame:
    rows: list[dict] = []

    for gene in genes:
        term = (
            f'{gene}[gene] AND ("pathogenic"[clinsig] OR "likely pathogenic"[clinsig]) '
            f'AND "p."[hgvs]'
        )
        s_params = {
            "db": "clinvar",
            "term": term,
            "retmax": max_per_gene,
            "retmode": "json",
        }
        s_resp = SESSION.get(ESEARCH_URL, params=s_params, timeout=30)
        s_resp.raise_for_status()
        id_list = s_resp.json().get("esearchresult", {}).get("idlist", [])

        time.sleep(0.4)
        if not id_list:
            continue

        batch_size = 50
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i : i + batch_size]
            sum_params = {
                "db": "clinvar",
                "id": ",".join(batch_ids),
                "retmode": "json",
            }
            sum_resp = SESSION.get(ESUMMARY_URL, params=sum_params, timeout=30)
            sum_resp.raise_for_status()
            result_map = sum_resp.json().get("result", {})

            for uid in batch_ids:
                doc = result_map.get(uid)
                if not isinstance(doc, dict):
                    continue

                title = doc.get("title") or doc.get("variation_set", [{}])[0].get("variation_name", "")
                if not title:
                    continue

                rows.append(
                    {
                        "gene_symbol": gene,
                        "variant_id": str(doc.get("uid", uid)),
                        "variant_title": str(title),
                        "clinical_significance": _extract_clinical_significance(doc),
                    }
                )

            time.sleep(0.4)

    return pd.DataFrame(rows).drop_duplicates(subset=["gene_symbol", "variant_id"]).reset_index(drop=True)


def main() -> None:
    print("Fetching UniProt proteins and domains...")
    proteins_df, domains_df = get_uniprot_data(GENES)
    proteins_df.to_csv(RAW_DIR / "proteins.csv", index=False)
    domains_df.to_csv(RAW_DIR / "domains.csv", index=False)

    print("Fetching ClinVar variants...")
    variants_df = get_clinvar_variants(GENES, max_per_gene=150)
    variants_df.to_csv(RAW_DIR / "variants.csv", index=False)

    print(f"Done! Proteins: {len(proteins_df)}, Domains: {len(domains_df)}, Variants: {len(variants_df)}")


if __name__ == "__main__":
    main()