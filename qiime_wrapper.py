import os
import subprocess
import glob
import pandas as pd

def make_manifest(reads=None, outdir=None, archive=None, log=None):
    """create manifest file for qiime2 importing"""
    paths = sorted(glob.glob(os.path.join(reads, "*")))
    manifest = os.path.join(outdir, "qiime_manifest.tsv")

    data = {"sample-id": [], "forward-absolute-filepath": [], "reverse-absolute-filepath": []}

    sample_num = 1
    for i in range(0, len(paths), 2):
        data["sample-id"].append("sample-" + str(sample_num))
        data["forward-absolute-filepath"].append(paths[i])
        data["reverse-absolute-filepath"].append(paths[i+1])
        sample_num += 1


    man_df = pd.DataFrame(data)
    man_df.to_csv(manifest, sep="\t", index=None)

    subprocess.run(
        ["qiime", "tools", "import",
         "--type", "SampleData[PairedEndSequencesWithQuality]",
         "--input-path", manifest,
         "--input-format", "PairedEndFastqManifestPhred33V2",
         "--output-path", archive],
         stdout=log,
         stderr=log
         )

def main():
    project_dir = os.getcwd()
    outdir = os.path.join(project_dir, "output")

    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # make directory to store output

    log = open(os.path.join(outdir, "wrapper.log"), "w") # open log

    reads = os.path.join(project_dir, "data", "subset_reads")
    qiime_archive = os.path.join(outdir, "subset_reads.qza")

    make_manifest(reads, outdir, qiime_archive, log)

    trimmed_reads = os.path.join(outdir, "trimmed_reads.qza")

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG")

    # subprocess.run(["qiime", "cutadapt", "trim-paired", "--i-demultiplexed-sequences", qiime_archive, "--p-front-f", primers[0], "--p-front-r", primers[1], "--o-trimmed-sequences", trimmed_reads])

if __name__ == "__main__":
    main()