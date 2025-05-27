import os
import sys
import subprocess
import glob
import argparse
import pandas as pd

threads = "12"

def get_args(args=None):
    """read in command line arguments"""
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    return parser.parse_args(args)

def import_reads(reads=None, outdir=None, log=None):
    """create manifest file for qiime2 and import reads"""
    paths = sorted(glob.glob(os.path.join(reads, "*")))
    manifest = os.path.join(outdir, "qiime_manifest.tsv")
    archive = os.path.join(outdir, "reads.qza")

    data = {"sample-id": [], "forward-absolute-filepath": [], "reverse-absolute-filepath": []}

    sample_num = 1
    for i in range(0, len(paths), 2):
        data["sample-id"].append("sample-" + str(sample_num))
        data["forward-absolute-filepath"].append(paths[i])
        data["reverse-absolute-filepath"].append(paths[i+1])
        sample_num += 1

    man_df = pd.DataFrame(data)
    man_df.to_csv(manifest, sep="\t", index=None)

    print("importing reads...")

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "SampleData[PairedEndSequencesWithQuality]",
        "--input-path", manifest,
        "--input-format", "PairedEndFastqManifestPhred33V2",
        "--output-path", archive
    ])
    
    print(f"done\n")
    
    return archive
    
def trim_reads(reads=None, outdir=None, primers=None, log=None):
    """trim reads using cutadapt"""
    trimmed_reads = os.path.join(outdir, "trimmed_reads.qza")

    print("trimming reads...")

    subprocess.run([
        "qiime", "cutadapt", "trim-paired",
        "--i-demultiplexed-sequences", reads,
        "--p-front-f", primers[0],
        "--p-front-r", primers[1],
        "--p-cores", threads,
        "--o-trimmed-sequences", trimmed_reads
    ])
    
    print(f"done\n")

    return trimmed_reads

def denoise_reads(trimmed_reads=None, outdir=None, log=None):
    """denoise reads using dada2"""
    asv_seqs = os.path.join(outdir, "asv-seqs.qza")

    print("running dada2...")

    subprocess.run([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", trimmed_reads,
        "--p-trunc-len-f", "95",
        "--p-trunc-len-r", "95",
        "--p-n-threads", threads,
        "--o-representative-sequences", asv_seqs,
        "--o-table", os.path.join(outdir, "feature-table.qza"),
        "--o-denoising-stats", os.path.join(outdir, "dada2-stats.qza")
    ])
    
    print(f"done\n")
    
    return asv_seqs
    
def run_vsearch(asv_seqs=None, ref_seqs=None, ref_taxa=None, outdir=None, log=None):
    out_taxa = os.path.join(outdir, "vsearch_taxa.qza")
    top_hits = os.path.join(outdir, "vsearch_top_hits.qza")

    print("running vsearch...")

    subprocess.run([
        "qiime", "feature-classifier", "classify-consensus-vsearch",
        "--i-query", asv_seqs,
        "--i-reference-reads", ref_seqs,
        "--i-reference-taxonomy", ref_taxa,
        "--p-perc-identity", "1.0",
        "--p-min-consensus", "0.94",
        "--p-threads", threads,
        "--o-classification", out_taxa,
        "--o-search-results", top_hits
    ])
    
    print(f"done\n")
    
    # return

def main():
    args = get_args(sys.argv[1:]) # get command line arguments

    project_dir = os.getcwd()
    outdir = os.path.join(project_dir, "output")

    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # make directory to store output

    log = open(os.path.join(outdir, "wrapper.log"), "w") # open log

    reads = os.path.join(project_dir, args.input) # path to reads from arguments

    qiime_archive = import_reads(reads, outdir, log)

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG") # MAYBE TAKE THIS AS INPUT IDK

    trimmed_reads = trim_reads(qiime_archive, outdir, primers, log)

    asv_seqs = denoise_reads(trimmed_reads, outdir, log)

    # asv_seqs = os.path.join(outdir, "asv-sequences-0.qza")

    # these are generated from the database file using the format script in tools/
    ref_seqs = os.path.join(project_dir, "data", "database", "seq_ref_for_qiime_vsearch.qza")
    ref_taxa = os.path.join(project_dir, "data", "database", "taxa_ref_for_qiime_vsearch.qza")

    run_vsearch(asv_seqs, ref_seqs, ref_taxa, outdir, log)

    log.close()

    print("fin")

if __name__ == "__main__":
    main()