# change database file from tsv to fasta for use in vsearch
import pandas as pd

tax_db = pd.read_csv("test/database/12S_Kankakee_db_1117.tsv", sep="\t")

with open("test/database/db_formatted_for_vsearch.fasta", "w") as f:
    for row in tax_db.itertuples(): # iterate through df, extract ids and sequences and write in fasta format
        id = ">" + row.id
        seq = row.Sequence
        f.write(f"{id}\n")
        f.write(f"{seq}\n")