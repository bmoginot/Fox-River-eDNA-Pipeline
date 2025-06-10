# Fox-River-eDNA-Pipeline
master's degree, here i come :^D

## About
This pipeline processes and analyzes short-read environmental DNA reads. The goal of my project specifically is to classify the taxa in the Fox River at 3 dams and determine if there is a significant difference in the taxa above and below each dam. 

## Workflow

### Reading Input into QIIME2
1. The path to input directory is read from the command line. This dir should contain demultiplexed illumina reads, with each file corresponding to one sample. 
2. Reads are then assembled into a *manifest file*. This is a tsv that related an arbitrary sample ID to each forward and reverse read pair. This manifest is used by QIIME2 to compress each each fastq file into a *archive file*. All output from this pipeline will be QIIME2 archive files, the contents of which can be viewed by using the bash *unzip* command.
3. Primers are then removed via Cutadapt.

### Denoising Reads
1. Trimmed reads are then passed to DADA2.
2. Reads are all truncated to a length of 95 bases at the filterAndTrim step. 
3. Reads are then denoised and pairs are merged to generate amplicon sequence variants (ASVs). 
4. Chimeras are removed.
5. Final output is (i) a table showing how many reads were removed at each step, (ii) a biom feature table, and (iii) a fasta file of each ASV detected by DADA2.

### Taxonomic Classification

#### Tools
This is a three-pronged approach to taxonomic classification.
1. VSEARCH looks for exact matches (100% identity) with 94% coverage.
2. A naive Bayes classifier is trained on the reference database and then used to classify sequences.
3. Nucleotide BLAST is a last resort to classify sequences that might be absent from the database (these sequences are likely to be contaminants such as human DNA).

#### Procedure
1. At each step, reads are compared against a pre-built database of 12S mitochondrial rRNA sequences. 
2. Sequences that are classified at the family level are set aside as "classified" sequences. The rest remain as "unclassified" sequences. Both are written out to a tsv.
3. The full ASV fasta from DADA2 is read in and implemented as a hash map. The feature IDs from the unclassified tsv are used to index this hash map to recover the sequences that are yet to be classified (since that information is not in the tsv). These sequences are complied in a fasta file and passed as input to the next tool.

## temp info
Here's a rough outline of the pipeline so far (for me)

1. input demultiplexed illumina reads
2. create qiime manifest based on filenames for reads, and use that manifest to create qiime archive for reads
3. trim reads with cutadapt
4. run dada2 and get fasta file (seqs are asvs)
5. pass fasta to vsearch (with input generated from database formatting script)
6. filter vserach output and pass to naive bayesian classifier

input files (qiime archives unless otherwise stated):
importing           manifest.tsv file with paths to input reads (created within the function call)
cutadapt            imported reads
dada2               trimmed reads
vsearch             asv sequences from dada2 (in fasta format)
                    fasta file with sequences and sequence ids
                    tsv file with sequence ids and taxonomic classification
naive bayes
    training        same as vsearch sans asv seqs
    classification  filtered asv seqs + classifier created in training

output files (output as qiime archives):
importing           directory of fastq.gz files
cutadapt            directory of fastq.gz files (sans primers)
dada2               fasta file; i believe each sequence corresponds to an asv but that info is not in this file
                    biom file containing a "feature table"; i'm not sure how to look at or interact with this
                    tsv file showing reads lost overall for each step        
vserach             tsv mapping feature id (based on dada2 fasta file names) to taxon
                    blast6 format tsv mapping feature id to sequence id (from kankakee database) + other stats
naive bayes
    training        tarred pickle file containing saved model
                    html file showing accuracy of model over the course of training (?)
                    tsv similar to the vsearch output (feature id is sequence name b/c its provided unlike with classification)
    classification  tsv similar to vsearch output (but for input sequences this time)