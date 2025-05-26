import os
import sys
import subprocess
import glob
import argparse
import pandas as pd

def get_args(args=None):
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    return parser.parse_args(args)

def make_manifest(reads=None, outdir=None, log=None):
    """create manifest file for qiime2 importing"""
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

    subprocess.run(
        ["qiime", "tools", "import",
         "--type", "SampleData[PairedEndSequencesWithQuality]",
         "--input-path", manifest,
         "--input-format", "PairedEndFastqManifestPhred33V2",
         "--output-path", archive],
         stdout=log,
         stderr=log
         )
    
    return archive
    
def trim_reads(reads=None, outdir=None, primers=None, log=None):
    trimmed_reads = os.path.join(outdir, "trimmed_reads.qza")

    subprocess.run(
        ["qiime", "cutadapt", "trim-paired",
         "--i-demultiplexed-sequences", reads,
         "--p-front-f", primers[0],
         "--p-front-r", primers[1],
         "--o-trimmed-sequences", trimmed_reads],
         stdout=log,
         stderr=log)

    return trimmed_reads

def main():
    project_dir = os.getcwd()
    outdir = os.path.join(project_dir, "output")

    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # make directory to store output

    log = open(os.path.join(outdir, "wrapper.log"), "w") # open log

    args = get_args(sys.argv[1:])
    reads = os.path.join(project_dir, args.input)

    qiime_archive = make_manifest(reads, outdir, log)

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG")

    trimmed_reads = trim_reads(qiime_archive, outdir, primers, log)

if __name__ == "__main__":
    main()