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

def run_dada2(reads=None, outdir=None, log=None, fasta_name=None):
    """
    denoises trimmed reads using dada2
    saves asv table as .tsv for downstream taxonomic classification
    """
    subprocess.run(
        ["Rscript", "src/denoise_reads.R", "-i", reads, "-o", outdir, "-f", fasta_name],
        stdout=log,
        stderr=log
        )
    
    return os.path.join(outdir, fasta_name)

def run_vsearch(outdir=None, query=None, seqs=None, taxa=None, log=None):
    """finds exact matches with vsearch"""
    out = os.path.join(outdir, "vserach_hits.tsv")
    tax_map = os.path.join(outdir, "vsearch_taxa.tsv")

    subprocess.run(
        ["vsearch", "--usearch_global", query,
         "--db", seqs,
         "--id", "1.0",
         "--query_cov", "0.94",
         "--strand", "plus",
         "--blast6out", out,
         "--maxaccepts", "0",
         "--maxhits", "1"],
         stdout=log,
         stderr=log
         )

    # map taxonomy using vsearch output and taxa info from database
    hits = pd.read_csv(out, sep="\t", header=None, names=[
        "query", "id", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ])
    taxonomy = pd.read_csv(taxa, sep="\t", header=None, names=["id", "taxon"])

    # Merge on reference ID
    merged = hits.merge(taxonomy, on="id", how="left")

    # Save to new file
    merged[["query", "taxon"]].to_csv(tax_map, sep="\t", index=False)

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

    fasta_name = "asv.fasta" # name of asv table to store for later

    asv_fasta = run_dada2(trimmed_reads_dir, outdir, log, fasta_name)

    seqs = os.path.join(project_dir, "data", "database", "vsearch_ref_seqs.fasta") # CHANGE THIS
    taxa = os.path.join(project_dir, "data", "database", "vsearch_ref_taxa.tsv") # CHANGE THIS
    run_vsearch(outdir, asv_fasta, seqs, taxa, log)

    log.close()

if __name__ == "__main__":
    main()