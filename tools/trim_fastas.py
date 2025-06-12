# remove unclassified seqs from dada2 fasta output (i <3 hashmaps)

import os
import glob

os.mkdir("tmp")
os.system("unzip -qd tmp output/dada2/asv-seqs.qza")
dada2_seqs = glob.glob("tmp/*/data/*")[0]

asv_map = {}

with open(dada2_seqs) as f: # read in dada2 asv seqs and construct hashmap with id: seq
    lines = f.readlines()
    for i in range(0, len(lines), 2):
        feat = lines[i][1:].strip()
        seq = lines[i+1].strip()
        asv_map[feat] = seq

os.system("rm -r tmp")

unclass_feats = []

with open(os.path.join("output", "bayes", "unassigned_bayes_seqs.fasta")) as f: # get all features left unclassified even after bayes
    lines = f.readlines()
    for i in range(0, len(lines), 2):
        feat = lines[i][1:].strip()
        unclass_feats.append(feat)

with open(os.path.join("phyloseq_input", "final_fasta.fasta"), "w") as f: # write out all seqs except unclassified ones
    for feat in asv_map.keys():
        if feat not in unclass_feats:
            f.write(f">{feat}\n")
            f.write(f"{asv_map[feat]}\n")