# Fox-River-eDNA-Pipeline
Here's a rough outline of the pipeline so far

1. input demultiplexed illumina reads
2. create qiime manifest based on filenames for reads, and use that manifest to create qiime archive for reads
3. trim reads with cutadapt
4. run dada2 and get fasta file (seqs are asvs)
5. pass fasta to vsearch
6. .

current output (output as qiime archives):
importing           directory of fastq.gz files
cutadapt            directory of fastq.gz files (sans primers)
dada2               fasta file; i believe each sequence corresponds to an asv but that info is not in this file
                    biom file containing a "feature table"; i'm not sure how to look at or interact with this
                    tsv file showing reads lost overall for each step        
vserach             tsv mapping feature id (based on dada2 fasta file names) to taxon
                    blast6 format tsv mapping feature id to sequence id (from kankakee database) + other stats
naive bayes         .

input files (qiime archives unless otherwise stated):
importing           manifest.tsv file with paths to input reads
cutadapt            imported reads
dada2               trimmed reads
vsearch             asv sequences from dada2 (in fasta format)
                    fasta file with sequences and sequence ids
                    tsv file with sequence ids and taxonomic classification
naive bayes         .