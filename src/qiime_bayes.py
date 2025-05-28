# train naive bayesian classifier and apply to dada2 output
import subprocess
import os
import pandas as pd

def unzip_qza(infile=None, outfile=None):
    """extract file from qiime archive"""
    tmp = os.path.join(os.getcwd(), "tmp")
    if not os.path.isdir("tmp"):
        os.system(f"mkdir {tmp}") # make tmp directory to store unzipped archive

    archive = os.path.join(tmp, "archive")
    os.system(f"unzip -q {infile} -d {archive}") # unzip qiime archive
    os.system(f"cp {archive}/*/data/* {outfile}") # remove file of interest from archive

    os.system(f"rm -r {archive}") # clean up temporary archive file

def parse_output(taxa_in=None):
    tmp = os.path.join(os.getcwd(), "tmp")
    if not os.path.isdir("tmp"):
        os.system(f"mkdir {tmp}") # create tmp directory to store taxa files

    unzipped_taxa = os.path.join(tmp, "taxa.tsv")
    unzip_qza(taxa_in, unzipped_taxa) # extract taxonomy tsv from qiime archive

    taxa_vsearch = pd.read_csv(unzipped_taxa, sep="\t") # read in vsearch output taxonomy.tsv

    unassigned = taxa_vsearch[
        taxa_vsearch["Taxon"] # take the column with taxonomic classification
        .str.split(";", expand=False) # split it into a list of taxonomic levels
        .apply(lambda x: len(x) < 5) # save the rows where there are fewer than 5 levels (less than family-level)
        ]
    
    unassigned.to_csv(os.path.join(tmp, "unassigned_vsearch_taxa.tsv"), sep="\t", index=False)
    
    retained = taxa_vsearch.drop( # drop all rows shared with unassigned (retain only family-level classifications)
        taxa_vsearch[
            taxa_vsearch["Feature ID"]
            .isin(unassigned["Feature ID"])
        ]
        .index
    )
    
    retained.to_csv(os.path.join(tmp, "retained_vsearch_taxa.tsv"), sep="\t", index=False)

    # os.system("rm -r tmp") # remove tmp directory

def nb_classifier():
    cwd = os.getcwd()
    train_seqs = os.path.join(cwd, "data", "database", "seq_ref_for_qiime_vsearch.qza")
    train_taxa = os.path.join(cwd, "data", "database", "taxa_ref_for_qiime_vsearch.qza")

    rescript_classifier = os.path.join(cwd, "test_out", "rescript_classifier")
    rescript_evaluation = os.path.join(cwd, "test_out", "rescript_evaluation")
    rescript_observed_taxonomy = os.path.join(cwd, "test_out", "rescript_observed_taxonomy")

    threads = "12"

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

    asv_seqs = os.path.join(cwd, "output", "asv-seqs.qza") # CHANGE
    nb_model = os.path.join(cwd, "test_out", "rescript_classifier.qza")
    nb_classification = os.path.join(cwd, "test_out", "nb_classification")

    print("classifying...")

    subprocess.run([
        "qiime", "feature-classifier", "classify-sklearn",
        "--i-reads", asv_seqs,
        "--i-classifier", nb_model,
        "--p-n-jobs", threads,
        "--o-classification", nb_classification
    ])

    print(f"done\n")

def main():
    vsearch_out = os.path.join(os.getcwd(), "output", "vsearch_taxa.qza")
    parse_output(vsearch_out)

if __name__ == "__main__":
    main()