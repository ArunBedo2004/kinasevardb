CREATE TABLE IF NOT EXISTS proteins (
    gene_symbol TEXT PRIMARY KEY,
    uniprot_id TEXT UNIQUE NOT NULL,
    protein_name TEXT NOT NULL,
    sequence_length INTEGER
);

CREATE TABLE IF NOT EXISTS domains (
    domain_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    uniprot_id TEXT NOT NULL,
    domain_name TEXT NOT NULL,
    start_pos INTEGER NOT NULL,
    end_pos INTEGER NOT NULL,
    FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol)
);

CREATE TABLE IF NOT EXISTS variants (
    variant_id TEXT PRIMARY KEY,
    gene_symbol TEXT NOT NULL,
    variant_title TEXT,
    clinical_significance TEXT,
    aa_position INTEGER,
    in_kinase_domain BOOLEAN,
    matched_domain TEXT,
    FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol)
);

CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_variants_domain ON variants(in_kinase_domain);
CREATE TABLE IF NOT EXISTS drugs (
    drug_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_symbol TEXT NOT NULL,
    chembl_molecule_id TEXT NOT NULL,
    drug_name TEXT NOT NULL,
    max_phase INTEGER,
    mechanism_of_action TEXT,
    FOREIGN KEY (gene_symbol) REFERENCES proteins(gene_symbol)
);

CREATE INDEX IF NOT EXISTS idx_drugs_gene ON drugs(gene_symbol);