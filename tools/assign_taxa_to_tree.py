import sys
import argparse
import os
import glob

def get_args(args):
	"""pares command line arguments"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", "--tree", help="tree file in newick format", required=True)
	parser.add_argument("-a", "--alignment", help="aligned sequences", required=True)
	parser.add_argument("-t", "--taxonomy", help="taxonomy file from qiime2", required=True)
	return parser.parse_args()

def extract_archive(arc, out):
	"""extract pertinent files from their qiime2 archive"""
	os.system(f"unzip -q {arc} -d tmp")
	target_file = glob.glob("tmp/*/data/*")[0]
	os.system(f"mv {target_file} {out}")
	os.system("rm -r tmp")
	return os.path.join(out, os.path.split(target_file)[-1])

def replace_taxonomy(tree, seqs, taxa):
	"""replace qiime ids in tree and alignment files with taxonomic assignments"""
	tax_hash = {} # hashmap to store taxonomy assignments

	with open(taxa) as f:
		next(f) # skip header
		for line in f.readlines(): # build hashmap where id = taxonomy
			fields = line.split("\t") # divide tsv into columns (qiime id, taxonomy, and some metric i forgot)
			qiime_id = fields[0] # qiime id assigned to ASV (maps back to sequence in alignment)
			qiime_taxa = fields[1].split(";")[-2:] # get the last two taxonomic levels (usually genus and species)
			tax_hash[qiime_id] = qiime_taxa # populate hashmap

	for qid, qtaxa in tax_hash.items():
		if qtaxa[-1].startswith("s__"): # if the assignment goes all the way to species level, the final assignment is Genus species
			final_taxa_assignment = " ".join([x[3:] for x in qtaxa])
		else: # else just write out the Family or Genus if that is the highest resolution
			final_taxa_assignment = qtaxa[-1][3:]

		os.system(f"sed -i 's/{qid}/{final_taxa_assignment}/g' {seqs}") # replace qiime ids with respective taxonomy in the alignment file
		os.system(f"sed -i 's/{qid}/{final_taxa_assignment}/g' {tree}") # do the same for the newick tree file

def main():
	args = get_args(sys.argv[1:])

	out = "biota_newick_tree_files"
	if not os.path.isdir(out):
		os.mkdir(out)

	tree = extract_archive(args.tree, out)
	seqs = extract_archive(args.alignment, out)
	taxa = extract_archive(args.taxonomy, out)

	replace_taxonomy(tree, seqs, taxa)

if __name__ == "__main__":
	main() 
