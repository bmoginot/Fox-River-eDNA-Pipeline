# Fox-River-eDNA-Pipeline
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

### Phyloseq
The script src/process_for_phyloseq.py collects biom feature table output by dada2 and metadata table. It then concatenates the taxonomy output from vsearch and the bayes classifier before importing it back into QIIME2. It then looks at the ASV sequences from dada2 and gets rid of any that were unclassified or contaminants before also importing it. It then generates a rooted phylogenetic tree using these sequences. 

