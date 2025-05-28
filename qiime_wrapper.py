import os
import sys
import subprocess
import glob
import argparse
import pandas as pd

def get_args(args=None):
    """read in command line arguments"""
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    parser.add_argument("-t", "--threads", help="# cpu cores to use for applicable functions")
    return parser.parse_args(args)

def import_reads(reads=None, outdir=None):
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
    
def trim_reads(reads=None, outdir=None, primers=None, threads=1):
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

def denoise_reads(trimmed_reads=None, outdir=None):
    """denoise reads using dada2"""
    asv_seqs = os.path.join(outdir, "asv-seqs.qza")

    print("running dada2...")

    subprocess.run([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", trimmed_reads,
        "--p-trunc-len-f", "95",
        "--p-trunc-len-r", "95",
        "--o-representative-sequences", asv_seqs,
        "--o-table", os.path.join(outdir, "feature-table.qza"),
        "--o-denoising-stats", os.path.join(outdir, "dada2-stats.qza")
    ])
    
    print(f"done\n")
    
    return asv_seqs
    
def run_vsearch(asv_seqs=None, ref_seqs=None, ref_taxa=None, outdir=None, threads=1):
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
    
    # return None

def nb_classifier(asv_seqs=None, train_seqs=None, train_taxa=None, outdir=None, threads=1):
    rescript_classifier = os.path.join(outdir, "rescript_classifier")
    rescript_evaluation = os.path.join(outdir, "rescript_evaluation")
    rescript_observed_taxonomy = os.path.join(outdir, "rescript_observed_taxonomy")

    print("training...")

    # fit model to 12S database
    subprocess.run([
        "qiime", "rescript", "evaluate-fit-classifier",
        "--i-sequences", train_seqs,
        "--i-taxonomy", train_taxa,
        "--p-n-jobs", threads,
        "--o-classifier", rescript_classifier,
        "--o-evaluation", rescript_evaluation,
        "--o-observed-taxonomy", rescript_observed_taxonomy
    ])

    print(f"done\n")

    nb_model = os.path.join(outdir, "rescript_classifier.qza")
    nb_classification = os.path.join(outdir, "nb_classification")

    print("classifying...")

    # classify sequences using model built above
    subprocess.run([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", asv_seqs,
        "--i-classifier", nb_model,
        "--p-n-jobs", threads,
        "--o-classification", nb_classification
    ])

    print(f"done\n")

    # return None

def main():
    args = get_args(sys.argv[1:]) # get command line arguments

    threads = str(args.threads) if args.threads else "1"
    print(threads)

    project_dir = os.getcwd()
    outdir = os.path.join(project_dir, "output")

    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # make directory to store output

    # log = open(os.path.join(outdir, "wrapper.log"), "w") # open logS

    reads = os.path.join(project_dir, args.input) # path to reads from arguments

    qiime_archive = import_reads(reads, outdir)

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG") # MAYBE TAKE THIS AS INPUT IDK

    trimmed_reads = trim_reads(qiime_archive, outdir, primers, threads)

    asv_seqs = denoise_reads(trimmed_reads, outdir)

    # these are generated from the database file using the format script in tools/
    ref_seqs = os.path.join(project_dir, "data", "database", "seq_ref_for_qiime_vsearch.qza")
    ref_taxa = os.path.join(project_dir, "data", "database", "taxa_ref_for_qiime_vsearch.qza")

    # taxonomic classification
    run_vsearch(asv_seqs, ref_seqs, ref_taxa, outdir, threads)
    nb_classifier(asv_seqs, ref_seqs, ref_taxa, outdir, threads)

    # log.close()

    print("fin")

if __name__ == "__main__":
    main()