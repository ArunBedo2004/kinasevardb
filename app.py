import sqlite3
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

DB_PATH = Path("kinasevardb.sqlite")

st.set_page_config(page_title="KinaseVarDB", page_icon="🧬", layout="wide")
st.title("🧬 KinaseVarDB — Cancer Kinase Variants & Inhibitors")
st.caption("Exploring UniProt kinase domains, ClinVar pathogenic variants, and ChEMBL targeted drugs")

if not DB_PATH.exists():
    st.error("kinasevardb.sqlite not found. Run the scripts in etl/ first.")
    st.stop()


@st.cache_data
def load_tables(db_path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    conn = sqlite3.connect(db_path)
    try:
        proteins = pd.read_sql_query("SELECT * FROM proteins", conn)
        domains = pd.read_sql_query("SELECT * FROM domains", conn)
        variants = pd.read_sql_query("SELECT * FROM variants", conn)
        drugs = pd.read_sql_query("SELECT * FROM drugs", conn)
    finally:
        conn.close()
    return proteins, domains, variants, drugs


proteins_df, domains_df, variants_df, drugs_df = load_tables(str(DB_PATH))

selected_gene = st.sidebar.selectbox(
    "Select Kinase Gene",
    options=sorted(proteins_df["gene_symbol"].unique()),
)
domain_only = st.sidebar.checkbox("Only show variants in kinase domain", value=False)

prot = proteins_df[proteins_df["gene_symbol"] == selected_gene].iloc[0]
g_domains = domains_df[domains_df["gene_symbol"] == selected_gene]
g_variants = variants_df[variants_df["gene_symbol"] == selected_gene].copy()
if domain_only:
    g_variants = g_variants[g_variants["in_kinase_domain"] == 1]
g_drugs = drugs_df[drugs_df["gene_symbol"] == selected_gene].sort_values(
    by=["max_phase", "drug_name"], ascending=[False, True]
)

st.subheader(f"{prot['gene_symbol']} — {prot['protein_name']}")
st.write(f"**UniProt ID:** `{prot['uniprot_id']}` | **Length:** `{prot['sequence_length']}` aa")

m1, m2, m3 = st.columns(3)
m1.metric("Domains", len(g_domains))
m2.metric("Variants Shown", len(g_variants))
m3.metric("Targeted Drugs (Phase ≥ 2)", len(g_drugs))

# Domain & Variant Map
st.markdown("### Domain & Variant Map")
seq_len = int(prot["sequence_length"]) if pd.notna(prot["sequence_length"]) else 1000

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=[1, seq_len],
        y=[0, 0],
        mode="lines",
        line=dict(color="#64748b", width=6),
        name="Protein Backbone",
        hoverinfo="skip",
    )
)

for _, drow in g_domains.iterrows():
    if pd.notna(drow["start_pos"]) and pd.notna(drow["end_pos"]):
        d_start = float(drow["start_pos"])
        d_end = float(drow["end_pos"])
        d_name = str(drow["domain_name"])
        is_kinase = "kinase" in d_name.lower()

        fill = "rgba(245, 158, 11, 0.25)" if is_kinase else "rgba(59, 130, 246, 0.18)"
        border = "rgba(245, 158, 11, 0.85)" if is_kinase else "rgba(59, 130, 246, 0.65)"

        fig.add_vrect(
            x0=d_start,
            x1=d_end,
            fillcolor=fill,
            line_width=1,
            line_color=border,
        )

        mid_x = (d_start + d_end) / 2.0
        fig.add_trace(
            go.Scatter(
                x=[mid_x],
                y=[0.0],
                mode="markers",
                marker=dict(size=14, opacity=0),
                name=d_name,
                showlegend=False,
                hovertemplate=(
                    f"<b>{d_name}</b><br>"
                    f"Range: {int(d_start)}–{int(d_end)} aa<extra></extra>"
                ),
            )
        )

mapped_variants = g_variants.dropna(subset=["aa_position"]).copy()
if not mapped_variants.empty:
    mapped_variants["aa_position"] = pd.to_numeric(mapped_variants["aa_position"], errors="coerce")
    mapped_variants = mapped_variants.dropna(subset=["aa_position"])

    in_dom = mapped_variants[mapped_variants["in_kinase_domain"] == 1]
    out_dom = mapped_variants[mapped_variants["in_kinase_domain"] != 1]

    if not out_dom.empty:
        fig.add_trace(
            go.Scatter(
                x=out_dom["aa_position"],
                y=[0.15] * len(out_dom),
                mode="markers",
                marker=dict(size=9, color="#94a3b8"),
                name="Outside Kinase Domain",
                text=out_dom["variant_title"],
                hovertemplate="<b>%{text}</b><br>AA Position: %{x}<extra></extra>",
            )
        )

    if not in_dom.empty:
        fig.add_trace(
            go.Scatter(
                x=in_dom["aa_position"],
                y=[0.15] * len(in_dom),
                mode="markers",
                marker=dict(size=11, color="#ef4444"),
                name="In Kinase Domain",
                text=in_dom["variant_title"],
                hovertemplate="<b>%{text}</b><br>AA Position: %{x}<extra></extra>",
            )
        )

fig.update_yaxes(visible=False, range=[-0.3, 0.55])
fig.update_xaxes(title_text="Amino Acid Position", range=[0, seq_len])
fig.update_layout(height=230, margin=dict(l=20, r=20, t=20, b=30))
st.plotly_chart(fig, width="stretch")

left_col, right_col = st.columns(2)

with left_col:
    st.markdown("### Protein & Kinase Domains")
    st.dataframe(g_domains[["domain_name", "start_pos", "end_pos"]], width="stretch", hide_index=True)

    st.markdown("### Variants")
    display_variants = g_variants.copy()

    def format_aa_pos(val) -> str:
        if pd.isna(val) or str(val).strip().lower() in ("", "none", "nan", "n/a"):
            return "N/A"
        num = pd.to_numeric(val, errors="coerce")
        if pd.notna(num):
            return str(int(num)) if float(num).is_integer() else str(num)
        return str(val)

    display_variants["aa_position"] = display_variants["aa_position"].apply(format_aa_pos).astype(str)
    st.dataframe(
        display_variants[["variant_id", "variant_title", "aa_position", "in_kinase_domain", "matched_domain"]],
        width="stretch",
        hide_index=True,
    )

with right_col:
    st.markdown("### Targeted Inhibitors")
    st.dataframe(
        g_drugs[["drug_name", "chembl_molecule_id", "max_phase", "mechanism_of_action"]],
        width="stretch",
        hide_index=True,
    )

# Cross-Gene Variant Overview
with st.expander("Cross-Gene Variant Summary (All Kinases)", expanded=False):
    gene_order = sorted(proteins_df["gene_symbol"].unique())
    summary_rows = []
    for gene in gene_order:
        sub = variants_df[variants_df["gene_symbol"] == gene]
        total_v = len(sub)
        in_kd = int((sub["in_kinase_domain"] == 1).sum()) if "in_kinase_domain" in sub.columns else 0
        summary_rows.append({
            "gene_symbol": gene,
            "total_variants": total_v,
            "in_kinase_domain": in_kd,
            "outside_or_unmapped": total_v - in_kd,
        })
    summary_df = pd.DataFrame(summary_rows)

    bar_fig = go.Figure()
    bar_fig.add_trace(
        go.Bar(
            x=summary_df["gene_symbol"],
            y=summary_df["in_kinase_domain"],
            name="In Kinase Domain",
            marker_color="#ef4444",
        )
    )
    bar_fig.add_trace(
        go.Bar(
            x=summary_df["gene_symbol"],
            y=summary_df["outside_or_unmapped"],
            name="Outside / Unmapped",
            marker_color="#64748b",
        )
    )
    bar_fig.update_layout(
        barmode="stack",
        height=280,
        xaxis_title="Kinase Gene",
        yaxis_title="Variant Count",
        margin=dict(l=20, r=20, t=25, b=30),
    )
    st.plotly_chart(bar_fig, width="stretch")
    st.dataframe(summary_df, width="stretch", hide_index=True)