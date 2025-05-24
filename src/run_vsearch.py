import os
import pandas as pd

cwd = os.getcwd()

query = os.path.join(cwd, "eDNA_pipeline_output", "asv.fasta")
seqs = os.path.join(cwd, "test", "database", "vsearch_ref_seqs.fasta")
out = os.path.join(cwd, "eDNA_pipeline_output", "vserach_hits.tsv")
taxa = os.path.join(cwd, "test", "database", "vsearch_ref_taxa.tsv")
tax_map = os.path.join(cwd, "eDNA_pipeline_output", "vsearch_taxa.tsv")

os.system(f"vsearch --usearch_global {query} \
          --db {seqs} \
          --id 1.0 \
          --query_cov 0.94 \
          --strand plus \
          --blast6out {out} \
          --maxaccepts 0 \
          --maxhits 1")

hits = pd.read_csv(out, sep="\t", header=None, names=[
    "query", "id", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
])
taxonomy = pd.read_csv(taxa, sep="\t", header=None, names=["id", "taxon"])

# Merge on reference ID
merged = hits.merge(taxonomy, on="id", how="left")

# Save to new file
merged[["query", "taxon"]].to_csv(tax_map, sep="\t", index=False)
