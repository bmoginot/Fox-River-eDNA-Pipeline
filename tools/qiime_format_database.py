# format database for use with vsearch in qiime2

import pandas as pd
import subprocess
import os

tax_db = pd.read_csv("data/database/12S_Kankakee_db_1117.tsv", sep="\t") # read in database

taxa_fasta = os.path.join(os.getcwd(), "data", "database", "kankakee_ref_seqs.fasta")
taxa_tsv = os.path.join(os.getcwd(), "data", "database", "kankakee_ref_taxa.tsv")

with open(taxa_fasta, "w") as f:
    for row in tax_db.itertuples(): # iterate through df, extract ids and sequences and write in fasta format
        id = ">" + row.id
        seq = row.Sequence
        f.write(f"{id}\n")
        f.write(f"{seq}\n")

# drop sequence, write taxa map to tsv
tax_db_noseq = tax_db.drop(labels="Sequence", axis=1)
tax_db_noseq.to_csv(taxa_tsv, sep="\t", index=False, header=["Feature ID", "Taxon"])

fasta_qza = os.path.join(os.getcwd(), "data", "database", "seq_ref_for_qiime_vsearch.qza")
ref_qza = os.path.join(os.getcwd(), "data", "database", "taxa_ref_for_qiime_vsearch.qza")

subprocess.run([
    "qiime", "tools", "import",
    "--type", "FeatureData[Sequence]",
    "--input-path", taxa_fasta,
    "--output-path", fasta_qza
    ])

subprocess.run([
    "qiime", "tools", "import",
    "--type", "FeatureData[Taxonomy]",
    "--input-format", "TSVTaxonomyFormat",
    "--input-path", taxa_tsv,
    "--output-path", ref_qza
    ])
