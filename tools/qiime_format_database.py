# format database for use with vsearch in qiime2

import pandas as pd

tax_db = pd.read_csv("data/database/12S_Kankakee_db_1117.tsv", sep="\t") # read in database

# drop sequence, write taxa map to tsv
tax_db_noseq = tax_db.drop(labels="id", axis=1)
tax_db_noseq.to_csv("data/database/qiime_ref_taxa_for_vsearch.tsv", sep="\t",
                    index=False, header=["Feature ID", "Taxon"])

