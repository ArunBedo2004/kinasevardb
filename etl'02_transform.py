from __future__ import annotations

from pathlib import Path
import re
import pandas as pd

RAW_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")

# Matches protein ranges like p.Glu746_Ala750del or p.Val600_Lys601insArg
AA_RANGE_RE = re.compile(r"p\.[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}(\d+)")
# Matches single-residue protein changes and frameshifts like p.Leu858Arg, p.Tyr827fs, p.Ile941fs*21
AA_SINGLE_RE = re.compile(r"p\.[A-Z][a-z]{2}(\d+)")


def parse_aa_positions(text: object) -> tuple[int | None, int | None]:
    """Extract (aa_start, aa_end) from an HGVS protein notation string."""
    if text is None or pd.isna(text):
        return None, None

    s = str(text)
    range_match = AA_RANGE_RE.search(s)
    if range_match:
        start_pos = int(range_match.group(1))
        end_pos = int(range_match.group(2))
        if start_pos > end_pos:
            start_pos, end_pos = end_pos, start_pos
        return start_pos, end_pos

    single_match = AA_SINGLE_RE.search(s)
    if single_match:
        pos = int(single_match.group(1))
        return pos, pos

    return None, None


def _pick_variant_text(row: pd.Series) -> str:
    """Combine candidate text columns so HGVS protein notations can be found."""
    parts: list[str] = []
    for col in ("variant_title", "hgvs_p", "protein_change", "name", "title"):
        if col in row.index and pd.notna(row[col]):
            parts.append(str(row[col]))
    return " ".join(parts)


def annotate_variants(
    variants_df: pd.DataFrame, domains_df: pd.DataFrame
) -> pd.DataFrame:
    df = variants_df.copy()

    # Build lookup of kinase domains by gene_symbol
    domains_by_gene: dict[str, list[tuple[str, int, int]]] = {}
    if not domains_df.empty and {"gene_symbol", "start_pos", "end_pos"}.issubset(
        domains_df.columns
    ):
        for _, drow in domains_df.iterrows():
            gene = str(drow["gene_symbol"]).strip()
            if pd.isna(drow["start_pos"]) or pd.isna(drow["end_pos"]):
                continue
            d_start = int(drow["start_pos"])
            d_end = int(drow["end_pos"])
            d_name = (
                str(drow["domain_name"])
                if "domain_name" in drow.index and pd.notna(drow["domain_name"])
                else "Kinase domain"
            )
            domains_by_gene.setdefault(gene, []).append((d_name, d_start, d_end))

    aa_positions: list[int | None] = []
    aa_end_positions: list[int | None] = []
    in_kinase_domain: list[int] = []
    matched_domain: list[str | None] = []

    for _, row in df.iterrows():
        raw_text = _pick_variant_text(row)
        start_pos, end_pos = parse_aa_positions(raw_text)

        # Keep existing aa_position if regex didn't find one
        if start_pos is None and "aa_position" in df.columns and pd.notna(row["aa_position"]):
            try:
                start_pos = int(float(row["aa_position"]))
                end_pos = start_pos
            except (ValueError, TypeError):
                start_pos, end_pos = None, None

        aa_positions.append(start_pos)
        aa_end_positions.append(end_pos)

        gene = str(row["gene_symbol"]).strip() if "gene_symbol" in df.columns and pd.notna(row["gene_symbol"]) else ""
        gene_domains = domains_by_gene.get(gene, [])

        found_domain: str | None = None
        if start_pos is not None and end_pos is not None:
            for d_name, d_start, d_end in gene_domains:
                if d_start <= end_pos and d_end >= start_pos:
                    found_domain = d_name
                    break

        in_kinase_domain.append(1 if found_domain is not None else 0)
        matched_domain.append(found_domain)

    df["aa_position"] = aa_positions
    df["aa_end_position"] = aa_end_positions
    df["in_kinase_domain"] = in_kinase_domain
    df["matched_domain"] = matched_domain
    return df


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    proteins_df = pd.read_csv(RAW_DIR / "proteins.csv")
    domains_df = pd.read_csv(RAW_DIR / "domains.csv")
    variants_df = pd.read_csv(RAW_DIR / "variants.csv", dtype={"variant_id": str})

    enriched_variants_df = annotate_variants(variants_df, domains_df)

    proteins_df.to_csv(PROCESSED_DIR / "clean_proteins.csv", index=False)
    domains_df.to_csv(PROCESSED_DIR / "clean_domains.csv", index=False)
    enriched_variants_df.to_csv(PROCESSED_DIR / "clean_variants.csv", index=False)

    parsed_count = int(enriched_variants_df["aa_position"].notna().sum())
    in_domain_count = int(enriched_variants_df["in_kinase_domain"].sum())

    print(f"Total variants: {len(enriched_variants_df)}")
    print(f"Variants with parsed AA position: {parsed_count}")
    print(f"Variants inside a kinase domain: {in_domain_count}")


if __name__ == "__main__":
    main()