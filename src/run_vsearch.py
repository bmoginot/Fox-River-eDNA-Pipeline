import os
import pandas as pd

cwd = os.getcwd()

query = os.path.join(cwd, "eDNA_pipeline_output", "asv.fasta")
seqs = os.path.join(cwd, "data", "database", "vsearch_ref_seqs.fasta")
out = os.path.join(cwd, "eDNA_pipeline_output", "vserach_hits.tsv")
taxa = os.path.join(cwd, "data", "database", "vsearch_ref_taxa.tsv")
tax_map = os.path.join(cwd, "eDNA_pipeline_output", "vsearch_taxa.tsv")

os.system(f"vsearch --usearch_global {query} \
          --db {seqs} \
          --id 1.0 \
          --query_cov 0.94 \
          --maxaccepts 10 \
          --maxrejects 32 \
          --blast6out {out} \
          --threads 6 \
          --top_hits_only")

# read in vsearch output
vsearch_hits = pd.read_csv(out, sep="\t", header=None, names=[
    "query", "id", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

# create dict that maps each asv to each exact match
lca_dict = {}
for row in vsearch_hits.itertuples():
    if row.query not in lca_dict:
        lca_dict[row.query] = [row.id]
    else:
        lca_dict[row.query].append(row.id)

ref = pd.read_csv(taxa, sep="\t") # read in tsv that maps id to taxa

for query in lca_dict.keys(): # iterate through asvs matched in vsearch
    consensus = ["k__Unclassified", "p__Unclassified", "c__Unclassified", "o__Unclassified", "f__Unclassified", "g__Unclassified", "s__Unclassified"] # default before finding lca
    matched_taxa = list(ref.loc[ref["id"].isin(lca_dict[query])]["Taxon"]) # get all rows in ref where ASV mapped to id
    if len(matched_taxa) > 1:
        split_taxa = [t.split(";") for t in matched_taxa] # divide string into individual taxa

        for i in range(len(split_taxa[0])): # compare each taxonomic level across each match for that asv
            level = set([x[i] for x in split_taxa])
            same = len(level) == 1
            if same:
                consensus[i] = next(iter(level)) # extract taxon level and replace at same level on consensus
    
    else:
        consensus = matched_taxa

    consensus = ";".join(consensus)
    print(query, consensus)