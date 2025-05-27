# train naive bayesian classifier and apply to dada2 output
import subprocess
import os

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

asv_seqs = os.path.join(cwd, "output", "asv-seqs.qza")
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