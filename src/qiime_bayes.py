import subprocess
import os

cwd = os.getcwd()
train_seqs = os.path.join(cwd, "data", "database", "seq_ref_for_qiime_vsearch.qza")
train_taxa = os.path.join(cwd, "data", "database", "taxa_ref_for_qiime_vsearch.qza")

rescript_classifier = os.path.join(cwd, "test_out", "rescript_classifier")
rescript_evaluation = os.path.join(cwd, "test_out", "rescript_evaluation")
rescript_observed_taxonomy = os.path.join(cwd, "test_out", "rescript_observed_taxonomy")

subprocess.run([
    "qiime", "rescript", "evaluate-fit-classifier",
    "--i-sequences", train_seqs,
    "--i-taxonomy", train_taxa,
    "--p-n-jobs", "12",
    "--o-classifier", rescript_classifier,
    "--o-evaluation", rescript_evaluation,
    "--o-observed-taxonomy", rescript_observed_taxonomy
    ])