import os
import pandas as pd
import glob
import subprocess

def format_metadata(dir):
    """get rid of guyana data and write metadata out to tsv"""
    print("formatting metadata...")

    reads = glob.glob("data/subset_reads/*")
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

def stitch_taxa(dir):
    """very simple command to create aggregate taxa file for phyloseq"""
    print ("stitching taxa files...")

    final_taxa = os.path.join(dir, "final_taxa.tsv")
    os.system(f"cp output/vsearch/retained_vsearch_taxa.tsv {final_taxa}") # duplicate vsearch output
    os.system(f"tail -n +2 output/bayes/retained_bayes_taxa.tsv >> {final_taxa}") # concatenate all but header from bayes output

    subprocess.run([
        "qiime", "tools", "import",
        "--type", "FeatureData[Taxonomy]",
        "--input-path", final_taxa,
        "--output-path", os.path.join(dir, "final_taxa")
    ])

    # os.remove(final_taxa)

    print(f"done\n")

def extract_qza(file, dir):
    """unzip qiime archive for formatting"""
    os.mkdir("tmp")
    os.system(f"unzip -qd tmp {file}")
    unzipped = glob.glob("tmp/*/data/*")[0]
    os.system(f"mv {unzipped} {dir}")
    os.system("rm -r tmp")

    return os.path.join(dir, os.path.split(unzipped)[-1])

def trim_fastas(dir):
    """grab fasta output from dada2 and remove sequences that were left unclassified after taxonomy steps"""
    print("trimming fasta file...")

    file = os.path.join("output", "dada2", "asv-seqs.qza")
    dada2_seqs = extract_qza(file, dir)

    asv_map = {}

    with open(dada2_seqs) as f: # read in dada2 asv seqs and construct hashmap with id: seq
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            feat = lines[i][1:].strip()
            seq = lines[i+1].strip()
            asv_map[feat] = seq

    # os.remove(dada2_seqs)

    unclass_feats = []

    with open(os.path.join("output", "bayes", "unassigned_bayes_seqs.fasta")) as f: # get all features left unclassified after bayes
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

    # os.remove(final_seqs)

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
    dir = "phyloseq_input"
    # if os.path.isdir(dir):
    #     os.system(f"rm -r {dir}")
    # os.mkdir(dir)

    os.system(f"cp output/dada2/feature-table.qza {dir}") # move feature table from dada2 output
    format_metadata(dir)
    # stitch_taxa(dir)
    # fasta = trim_fastas(dir)
    # align_to_tree(fasta, dir)

    # os.system(f"cp -r {dir} /mnt/c/Users/bmogi/OneDrive/Documents/UniDocs/MSThesis/") # export for analysis in R

if __name__ == "__main__":
    main()