# get rid of kankakee data and write metadata out to tsv to prep for dada2 run
import pandas as pd

metadata = pd.read_csv("kankakeerun_metadata_forqiime.txt", sep="\t")
metadata.head()

only_kankakee_data = metadata[metadata["Study"]=="Kankakee"]
only_kankakee_data.to_csv("only_kankakee_metadata.tsv", sep="\t")