import os
import sys
import subprocess
import glob
import argparse
import pandas as pd
import time

def get_args(args):
    """read in command line arguments"""
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    parser.add_argument("-t", "--threads", help="# cpu cores to use for applicable functions")
    return parser.parse_args(args)

def import_reads(reads, outdir):
    """create manifest file for qiime2 and import reads"""

    print("making manifest...")

    paths = sorted(glob.glob(os.path.join(reads, "*")))
    manifest = os.path.join(outdir, "qiime_manifest.tsv")
    archive = os.path.join(outdir, "reads.qza")

    data = {"sample-id": [], "forward-absolute-filepath": [], "reverse-absolute-filepath": []}

    for i in range(0, len(paths), 2):
        samid = paths[i].split("/")[-1].split("_")[0] # get the sample id from the path (corresponds to metadata)
        data["sample-id"].append(samid)
        data["forward-absolute-filepath"].append(paths[i])
        data["reverse-absolute-filepath"].append(paths[i+1])

    man_df = pd.DataFrame(data)
    man_df.to_csv(manifest, sep="\t", index=False)

    print(f"done\n")

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
    
def trim_reads(reads, outdir, primers, threads):
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

def denoise_reads(trimmed_reads, outdir):
    """denoise reads using dada2"""
    outdir = os.path.join(outdir, "dada2")
    os.mkdir(outdir)

    asv_seqs = os.path.join(outdir, "asv-seqs.qza") # goes into taxonomic classification
    feat_table = os.path.join(outdir, "feature-table.qza") # used for phyloseq

    print("running dada2...")

    subprocess.run([
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", trimmed_reads,
        "--p-trunc-len-f", "95",
        "--p-trunc-len-r", "95",
        "--o-representative-sequences", asv_seqs,
        "--o-table", feat_table,
        "--o-denoising-stats", os.path.join(outdir, "dada2-stats.qza")
    ])
    
    print(f"done\n")
    
    return asv_seqs, feat_table

def extract_qza(file, dir):
    """unzip qiime archive for formatting"""
    os.mkdir("tmp")
    os.system(f"unzip -qd tmp {file}")
    unzipped = glob.glob("tmp/*/data/*")[0]
    os.system(f"mv {unzipped} {dir}")
    os.system("rm -r tmp")

    return os.path.join(dir, os.path.split(unzipped)[-1])

def parse_output(taxa_in, taxa_out, outdir):
    """parse output from taxa classification and separate into retained and unassigned tsvs"""
    unzipped_taxa = extract_qza(taxa_in, outdir) # extract taxonomy tsv from qiime archive

    taxa_vsearch = pd.read_csv(unzipped_taxa, sep="\t") # read in vsearch output taxonomy.tsv
    os.remove(unzipped_taxa)

    outdir = os.path.join(outdir, taxa_out)

    print(f"parsing {taxa_out} taxonomy...")

    unassigned = taxa_vsearch[
        taxa_vsearch["Taxon"] # take the column with taxonomic classification
        .str.split(";", expand=False) # split it into a list of taxonomic levels
        .apply(lambda x: len(x) < 5) # save the rows where there are fewer than 5 levels (less than family-level)
        ]
    
    unassigned_out = os.path.join(outdir, "unassigned_" + taxa_out + "_taxa.tsv") # this goes into the map_seqs function to figure out which exact seqs were left unassigned --> bayes classifier
    unassigned.to_csv(os.path.join(unassigned_out), sep="\t", index=False)
    
    retained = taxa_vsearch.drop( # drop all rows shared with unassigned (retain only family-level classifications)
        taxa_vsearch[
            taxa_vsearch["Feature ID"]
            .isin(unassigned["Feature ID"])
        ]
        .index
    )
    
    retained_out = os.path.join(outdir, "retained_" + taxa_out + "_taxa.tsv") # this is used to prep taxonomy file for phyloseq
    retained.to_csv(retained_out, sep="\t", index=False)

    print(f"done\n")

    return unassigned_out, retained_out

def map_seqs(asv_seqs, unassigned, taxa_out, outdir):
    """recover fasta file after filtering out classified sequences"""
    unzipped_asvs = extract_qza(asv_seqs, outdir)

    print("mapping seqs...")

    asv_map = {}

    with open(unzipped_asvs) as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            feat = lines[i][1:].strip()
            seq = lines[i+1].strip()
            asv_map[feat] = seq
    
    os.remove(unzipped_asvs)

    features_index = list(pd.read_csv(unassigned, sep="\t")["Feature ID"])

    outdir = os.path.join(outdir, taxa_out)
    
    unassigned_fasta = os.path.join(outdir, "unassigned_" + taxa_out + "_seqs.fasta") # non-archive version of below which i pass to phyloseq processing part to trim taxa that are still unassigned after bayesian classification
    with open(unassigned_fasta, "w") as f:
        for feat in features_index:
            f.write(f">{feat}\n")
            f.write(f"{asv_map[feat]}\n")

    asv_out = os.path.join(outdir, "unassigned_" + taxa_out + "_seqs.qza") # unassigned seqs that go to the next classifier 

    subprocess.run([
        "qiime", "tools", "import",
        "--input-path", unassigned_fasta,
        "--output-path", asv_out,
        "--type", "FeatureData[Sequence]"
    ])

    print(f"done\n")

    return asv_out, unassigned_fasta

def run_vsearch(asv_seqs, ref_seqs, ref_taxa, outdir, threads):
    """run vsearch to classify taxa based on reference database"""
    outdir = os.path.join(outdir, "vsearch")
    os.mkdir(outdir)
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
    
    return out_taxa

def nb_classifier(asv_seqs, train_seqs, train_taxa, outdir, threads):
    """train naive bayesian model on reference database and classify sequences missed by vsearch"""
    outdir = os.path.join(outdir, "bayes")
    os.mkdir(outdir)
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
    nb_classification = os.path.join(outdir, "nb_classification.qza")

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

    return nb_classification

def format_metadata(dir, reads):
    """get rid of guyana data and write metadata out to tsv"""
    print("formatting metadata...")

    reads = glob.glob(os.path.join(reads, "*"))
    samples = [x.split("/")[-1].split("_")[0] for x in reads] # get the sample name (first field of illumina designation) for only the reads i'm using in my test

    metadata = pd.read_csv(os.path.join("data", "metadata", "kankakeerun_metadata_forqiime.txt"), sep="\t")

    subset_data = metadata[metadata["sampleID"].isin(samples)] # subset for only Kankakee data

    # metadata must be qiime2-compliant
    outfile = os.path.join(dir, "subset_kankakee_metadata.tsv")
    subset_data.to_csv(outfile, sep="\t", index=False)
    os.system(f"sed -i '1s/^sampleID/#SampleID/' {outfile}") # change first column header
    os.system(f"sed -i '1a \
              #q2:types\tnumeric\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical' \
              {outfile}") # add type annotations for each column

    print(f"done\n")

def stitch_taxa(dir, vsearch, bayes):
    """very simple command to create aggregate taxa file for phyloseq"""
    print ("stitching taxa files...")

    final_taxa = os.path.join(dir, "final_taxa.tsv")
    os.system(f"cp {vsearch} {final_taxa}") # duplicate vsearch output
    os.system(f"tail -n +2 {bayes} >> {final_taxa}") # concatenate all but header from bayes output

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-path", final_taxa,
        "--output-path", os.path.join(dir, "final_taxa")
    ])

    os.remove(final_taxa)

    print(f"done\n")

def trim_fastas(dir, asv_seqs, unassigned_bayes):
    """grab fasta output from dada2 and remove sequences that were left unclassified after taxonomy steps"""
    print("trimming fasta file...")

    file = asv_seqs
    dada2_seqs = extract_qza(file, dir)

    asv_map = {}

    with open(dada2_seqs) as f: # read in dada2 asv seqs and construct hashmap with id: seq
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            feat = lines[i][1:].strip()
            seq = lines[i+1].strip()
            asv_map[feat] = seq

    os.remove(dada2_seqs)

    unclass_feats = []

    with open(unassigned_bayes) as f: # get all features left unclassified after bayes
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            feat = lines[i][1:].strip()
            unclass_feats.append(feat)

    final_seqs = os.path.join(dir, "final_seqs.fasta")

    with open(final_seqs, "w") as f: # write out all seqs except unclassified ones
        for feat in asv_map.keys():
            if feat not in unclass_feats:
                f.write(f">{feat}\n")
                f.write(f"{asv_map[feat]}\n")

    achive_out = os.path.join(dir, "final_seqs.qza")

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Sequence]",
        "--input-path", final_seqs,
        "--output-path", achive_out
    ])

    os.remove(final_seqs)

    print(f"done\n")

    return achive_out

def align_to_tree(file, dir):
    """run mafft and fasttree in qiime to generate rooted tree for phyloseq"""
    print("aligning to tree...")

    subprocess.run([
          "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
          "--i-sequences", file,
          "--p-n-threads", "12",
          "--o-alignment", os.path.join(dir, "unmaksed_alignment"),
          "--o-masked-alignment", os.path.join(dir, "masked_alignment"),
          "--o-tree", os.path.join(dir, "unrooted_tree"),
          "--o-rooted-tree", os.path.join(dir, "rooted_tree")
    ])

    print(f"done\n")

def main():
    start = time.time()

    args = get_args(sys.argv[1:]) # get command line arguments

    threads = str(args.threads) if args.threads else "1"

    project_dir = os.getcwd()
    outdir = os.path.join(project_dir, "output") # TRY: os.path.abspath("outdir")

    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # make directory to store output

    # log = open(os.path.join(outdir, "wrapper.log"), "w") # open log

    reads = os.path.join(project_dir, args.input) # path to reads from arguments

    qiime_archive = import_reads(reads, outdir)

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG") # MAYBE TAKE THIS AS INPUT IDK

    trimmed_reads = trim_reads(qiime_archive, outdir, primers, threads)

    asv_seqs, feat_table = denoise_reads(trimmed_reads, outdir)

    # these are generated from the database file using the format script in tools/
    ref_seqs = os.path.join(project_dir, "data", "database", "seq_ref_for_qiime_vsearch.qza")
    ref_taxa = os.path.join(project_dir, "data", "database", "taxa_ref_for_qiime_vsearch.qza")

    # taxonomic classification
    vsearch_out = run_vsearch(asv_seqs, ref_seqs, ref_taxa, outdir, threads)
    unassigned_vsearch_taxa, retained_vsearch_taxa = parse_output(vsearch_out, "vsearch", outdir)

    unassigned_vsearch_seqs_archive, unassigned_vsearch_fasta = map_seqs(asv_seqs, unassigned_vsearch_taxa, "vsearch", outdir)

    bayes_out = nb_classifier(unassigned_vsearch_seqs_archive, ref_seqs, ref_taxa, outdir, threads)
    unassigned_bayes_taxa, retained_bayes_taxa = parse_output(bayes_out, "bayes", outdir)

    unassigned_bayes_seqs_archive, unassigned_bayes_fasta = map_seqs(asv_seqs, unassigned_bayes_taxa, "bayes", outdir)

    # prep files for phyloseq
    dir = os.path.join(outdir, "phyloseq_input")
    os.mkdir(dir)

    os.system(f"cp {feat_table} {dir}") # move feature table from dada2 output

    format_metadata(dir, reads)

    stitch_taxa(dir, retained_vsearch_taxa, retained_bayes_taxa)

    fasta = trim_fastas(dir, asv_seqs, unassigned_bayes_fasta)
    align_to_tree(fasta, dir)

    os.system(f"cp -r {dir} /mnt/c/Users/bmogi/OneDrive/Documents/UniDocs/MSThesis/") # export for analysis in R

    end = time.time()

    print(f"took {end - start} seconds\n")

    # CHECK OUT: classify-consensus-blast

    # log.close()

    print("fin")

if __name__ == "__main__":
    main()