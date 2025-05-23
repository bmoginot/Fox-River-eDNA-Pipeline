import os
import sys
import argparse
import glob
import subprocess

def get_args(args=None):
    parser = argparse.ArgumentParser(description="run eDNA pipeline")
    parser.add_argument("-i", "--input", help="directory containing reads", required=True)
    return parser.parse_args(args)

def trim(reads=None, primers=None, log=None):
    """
    trim reads
    first, creates directory to store trimmed reads. then, iterates through input directory and pulls out each read pair.
    then, puts each read pair through cutadapt with input primers and error rate. cut adapt trims the primers from the reads
    and writes the new fq files to the newly-created directory.
    """
    os.mkdir("trimmed_reads")
    os.chdir("trimmed_reads")
    path = reads + "*"
    files = sorted(glob.glob(path))

    erate = "0.25" # needs to be a string for subprocess.run

    for i in range(0, len(files), 2): # REMOVE HARCODING LATER; LET THEM PICK PRIMERS AND ERROR RATE
        fread = files[i]
        rread = files[i+1]
        fout = fread.split("/")[-1].split(".")[0] + "-trimmed.fastq"
        rout = rread.split("/")[-1].split(".")[0] + "-trimmed.fastq"
        subprocess.run(
            ["cutadapt", "-g", primers[0], "-G", primers[1], "-e", erate, "-o", fout, "-p", rout, fread, rread],
            stdout=log # write output to log
        )

def run_dada2(reads=None):
    pass

def main():
    args = get_args(sys.argv[1:])
    cwd = os.getcwd()
    reads_dir = os.path.join(cwd, args.input) # store path to input data directory

    outdir = "eDNA_pipeline_output"
    if os.path.isdir(outdir):
        os.system(f"rm -r {outdir}")
    os.mkdir(outdir) # create output directory and move into it
    os.chdir(outdir)

    log = open("wrapper.log", "w") # open log

    primers = ("ACTGGGATTAGATACCCC", "TAGAACAGGCTCCTCTAG")
    trim(reads_dir, primers, log)

    # run_dada2(reads_dir)

    log.close()

if __name__ == "__main__":
    main()