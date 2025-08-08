import glob
import sys
import os

def extract_qza(arc, out):
	os.system(f"unzip {arc} -d tmp")

	if not os.path.isdir(out):
		os.mkdir(out)

	os.system(f"mv tmp/*/data/* {out}")
	os.system("rm -r tmp")

def main():
	archive = sys.argv[1]
	outdir = sys.argv[2]
	extract_qza(archive, outdir)

if __name__ == "__main__":
	main()
