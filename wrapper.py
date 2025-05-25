import os
import sys
import argparse
import glob
import subprocess
import pandas as pd

def get_args(args=None):
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    return parser.parse_args(args)

def trim(reads=None, primers=None, outdir=None, log=None):
    """
    trim reads
    first, creates directory to store trimmed reads. then, iterates through input directory and pulls out each read pair.
    then, puts each read pair through cutadapt with input primers and error rate. cut adapt trims the primers from the reads
    and writes the new fq files to the newly-created directory.
    """
    trimmed_reads_dir = os.path.join(outdir, "trimmed_reads")
    os.mkdir(trimmed_reads_dir)
    path = reads + "*"
    files = sorted(glob.glob(path))

    erate = "0.25" # needs to be a string for subprocess.run

    for i in range(0, len(files), 2): # REMOVE HARCODING LATER; LET THEM PICK PRIMERS AND ERROR RATE
        fread = files[i]
        rread = files[i+1]
        fout = os.path.join(trimmed_reads_dir, fread.split("/")[-1].split(".")[0] + "-trimmed.fastq.gz")
        rout = os.path.join(trimmed_reads_dir, rread.split("/")[-1].split(".")[0] + "-trimmed.fastq.gz")
        subprocess.run(
            ["cutadapt",
             "-g", primers[0],
             "-G", primers[1],
             "-e", erate,
             "-o", fout,
             "-p", rout,
             fread, rread],
            stdout=log, # write output to log
            stderr=log # and error
            )

    return trimmed_reads_dir

def run_dada2(reads=None, outdir=None, log=None):
    """
    denoises trimmed reads using dada2
    saves asv table as .tsv for downstream taxonomic classification
    """
    subprocess.run(
        ["Rscript", "src/denoise_reads.R", "-i", reads, "-o", outdir],
        stdout=log,
        stderr=log
        )

def run_vsearch(outdir=None, query=None, seqs=None, taxa=None, log=None):
    """finds exact matches with vsearch and computes lca"""
    out = os.path.join(outdir, "vserach_hits.tsv")
    tax_map = os.path.join(outdir, "vsearch_lca.tsv")

    subprocess.run(
        ["vsearch", "--usearch_global", query,
         "--db", seqs,
         "--id", "1.0",
         "--query_cov", "0.94",
         "--maxaccepts", "10",
         "--maxrejects", "32",
         "--blast6out", out,
         "--threads", "6",
         "--top_hits_only"],
         stdout=log,
         stderr=log
         )

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

    final_taxa_map = {"query": [], "taxon": []} # data for final tsv output

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
            consensus = matched_taxa # skip above if there was only one match

        final_taxa_map["query"].append(query) # asv number
        final_taxa_map["taxon"].append(";".join(consensus)) # lca taxonomy

    lca_df = pd.DataFrame(data = final_taxa_map, columns=["query", "taxon"])
    lca_df.to_csv(tax_map, sep="\t")

def main():
    args = get_args(sys.argv[1:])
    project_dir = os.getcwd()
    reads_dir = os.path.join(project_dir, args.input) # store path to input data directory

    outdir = os.path.join(project_dir, "eDNA_pipeline_output")
    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # create output directory and move into it

    log = open(os.path.join(outdir, "wrapper.log"), "w") # open log

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG")
    trimmed_reads_dir = trim(reads_dir, primers, outdir, log)

    asv_fasta = run_dada2(trimmed_reads_dir, outdir, log)

    seqs = os.path.join(project_dir, "data", "database", "vsearch_ref_seqs.fasta") # CHANGE THIS
    taxa = os.path.join(project_dir, "data", "database", "vsearch_ref_taxa.tsv") # CHANGE THIS
    run_vsearch(outdir, asv_fasta, seqs, taxa, log)

    log.close()

if __name__ == "__main__":
    main()