# generate manifest file to import data into qiime2
import os
import glob
import pandas as pd

subset_reads = os.path.join(os.getcwd(), "data", "subset_reads", "*")
manifest_file = os.path.join(os.getcwd(), "output", "kankakee_manifest.tsv")

paths = sorted(glob.glob(subset_reads))

data = {"sample-id": [], "forward_absolute_path": [], "reverse_absolute_path": []}

sample_num = 1
for i in range(0, len(paths), 2):
    data["sample-id"].append("sample-" + str(sample_num))
    data["forward_absolute_path"].append(paths[i])
    data["reverse_absolute_path"].append(paths[i+1])
    sample_num += 1

man_df = pd.DataFrame(data)
man_df.to_csv(manifest_file, sep="\t")