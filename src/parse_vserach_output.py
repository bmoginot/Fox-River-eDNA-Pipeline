# parse taxonomy.tsv output by vsearch
# filter out at least family-level classifications and save for later
# retain rest and pass only those seqs to the naive bayes classifier
# do the same for that output

import os
import pandas as pd

def unzip_qza(infile=None, outfile=None):
    """extract file from qiime archive"""
    os.system(f"unzip -q {infile} -d tmp")
    os.system(f"mv tmp/*/data/* {outfile}")
    os.system("rm -r tmp")

def parse_vsearch_output():
    cwd = os.getcwd()
    vsearch_out = os.path.join(cwd, "output", "vsearch_taxa.qza")
    unzipped_vsearch = os.path.join(cwd, "vsearch_taxa.tsv")

    unzip_qza(vsearch_out, unzipped_vsearch) # extract taxonomy tsv from qiime archive

    taxa_vsearch = pd.read_csv(unzipped_vsearch, sep="\t")
    print(taxa_vsearch.head())

parse_vsearch_output()