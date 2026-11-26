from __future__ import annotations

from pathlib import Path
import sqlite3
import pandas as pd

PROCESSED_DIR = Path("data/processed")
DB_PATH = Path("kinasevardb.sqlite")
SCHEMA_PATH = Path("schema.sql")


def _table_columns(conn: sqlite3.Connection, table_name: str) -> list[str]:
    rows = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
    return [row[1] for row in rows]


def _append_matching_columns(
    df: pd.DataFrame,
    conn: sqlite3.Connection,
    table_name: str,
    dedupe_col: str | None = None,
) -> None:
    cols = _table_columns(conn, table_name)
    if not cols:
        raise RuntimeError(f"Target SQLite table '{table_name}' does not exist.")

    filtered = df[[c for c in cols if c in df.columns]].copy()

    if "variant_id" in filtered.columns:
        filtered["variant_id"] = filtered["variant_id"].astype(str).str.strip()
        filtered = filtered[filtered["variant_id"] != ""]

    if dedupe_col and dedupe_col in filtered.columns:
        filtered = filtered.drop_duplicates(subset=[dedupe_col], keep="first")

    filtered.to_sql(table_name, conn, if_exists="append", index=False)


def main() -> None:
    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA foreign_keys = OFF;")
        if SCHEMA_PATH.exists():
            schema_sql = SCHEMA_PATH.read_text(encoding="utf-8")
            conn.executescript(schema_sql)

        conn.execute("DELETE FROM variants")
        conn.execute("DELETE FROM domains")
        conn.execute("DELETE FROM proteins")

        proteins_df = pd.read_csv(PROCESSED_DIR / "clean_proteins.csv")
        domains_df = pd.read_csv(PROCESSED_DIR / "clean_domains.csv")
        variants_df = pd.read_csv(PROCESSED_DIR / "clean_variants.csv", dtype={"variant_id": str})

        _append_matching_columns(proteins_df, conn, "proteins", dedupe_col="gene_symbol")
        _append_matching_columns(domains_df, conn, "domains")
        _append_matching_columns(variants_df, conn, "variants", dedupe_col="variant_id")

        conn.execute("PRAGMA foreign_keys = ON;")
        conn.commit()

        query = """
        SELECT
            p.gene_symbol,
            p.protein_name,
            COUNT(v.variant_id) AS total_variants,
            SUM(CASE WHEN v.in_kinase_domain THEN 1 ELSE 0 END) AS kinase_domain_variants
        FROM proteins p
        LEFT JOIN variants v ON p.gene_symbol = v.gene_symbol
        GROUP BY p.gene_symbol, p.protein_name
        ORDER BY kinase_domain_variants DESC
        """
        summary_df = pd.read_sql_query(query, conn)
        print("Loaded into kinasevardb.sqlite! Summary by gene:")
        print(summary_df.to_string(index=False))
    finally:
        conn.close()


if __name__ == "__main__":
    main()