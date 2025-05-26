# format database for use with vsearch in qiime2

import pandas as pd
import subprocess
import os

tax_db = pd.read_csv("data/database/12S_Kankakee_db_1117.tsv", sep="\t") # read in database

taxa_tsv = os.path.join(os.getcwd(), "data", "database", "kankakee_ref.tsv")

# drop sequence, write taxa map to tsv
tax_db_noseq = tax_db.drop(labels="id", axis=1)
tax_db_noseq.to_csv(taxa_tsv, sep="\t", index=False, header=["Feature ID", "Taxon"])
                    
ref_qza = os.path.join(os.getcwd(), "data", "database", "ref_for_qiime_vsearch.qza")

subprocess.run([
    "qiime", "tools", "import",
    "--type", "FeatureData[Taxonomy]",
    "--input-format", "TSVTaxonomyFormat",
    "--input-path", taxa_tsv,
    "--output-path", ref_qza])
