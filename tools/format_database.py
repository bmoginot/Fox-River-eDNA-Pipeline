# change database file from tsv to fasta for use in vsearch
import pandas as pd

tax_db = pd.read_csv("test/database/12S_Kankakee_db_1117.tsv", sep="\t")

with open("test/database/vsearch_ref_seqs.fasta", "w") as f:
    for row in tax_db.itertuples(): # iterate through df, extract ids and sequences and write in fasta format
        id = ">" + row.id
        seq = row.Sequence
        f.write(f"{id}\n")
        f.write(f"{seq}\n")

# drop sequence, write taxa map to tsv
tax_db_noseq = tax_db.drop(labels="Sequence", axis=1)
tax_db_noseq.to_csv("test/database/vsearch_ref_taxa.tsv", sep="\t", index=False)