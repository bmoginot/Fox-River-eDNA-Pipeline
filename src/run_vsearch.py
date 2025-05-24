import os

cwd = os.getcwd()

query = os.path.join(cwd, "eDNA_pipeline_output", "asv.fasta")
ref = os.path.join(cwd, "test", "database", "vsearch_ref_seqs.fasta")
out = "hits.tsv"

os.system(f"vsearch --usearch_global {query} \
          --db {ref} \
          --id 1.0 \
          --query_cov 0.94 \
          --strand plus \
          --blast6out {out} \
          --maxaccepts 0 \
          --maxhits 1")