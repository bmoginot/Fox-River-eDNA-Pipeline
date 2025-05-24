# affirmation: I AM GOING TO TAKE VERY GOOD, DETAILED NOTES

alright lets get into this shit

## 05/19/25
i have read through most of dr. picq's kankakee paper and i think i have a pretty good idea of the methods. in the end i think i'm going to use qiime but for the time being i'm not gonna touch that shit.

goal: get the dada2 pipeline up and running on a subset of kankakee data

plan of attack:
download dada2 code from posit *
download kankakee data *
subset above *
format and clean metadata *
put it all into dada2

what i did:
logged into posit cloud, copied code to run dada2 that i used in my final project for metagenomics
accessed onedrive folder for gyuana-kankakee data and downloaded all 3 replicates for the first 4 sites, all of which are at rock creek
each of these replicates has forward and reverse reads, for a total of 24 fq files
each fq file has somewhere in the ballpark of 50-60k reads
i subset each of them down to 10k reads (40k lines) using a custom bash script (subset.sh)
subset metadata to get rid of guyana info and wrote it to a tsv

plan for tomorrow:
run dada2 tutorial to familiarize myself with the pipeline again
put subset reads through pipeline
    maybe get some more at different sites
    make sure metadata works

## 05/21/25
it was rainy so i didn't get into the lab yesterday

plan for today:
run through dada2 pipeline *
subset kankakee reads so they are a similar size to the mothur test data *
get cutadapt running and put into pipeline
get dada2 script running with a function call

what i did:
dada2 downloaded (install via bioconductor in R); mothur testdata downloaded
ran through dada2 tutorial and took notes in the dada2_tutorial.R script
average size of mothur reads is about equal to size of already subset kankakee data so no change there
new test data is reps 1-3 for samples 1-3 and 5-7. former are from rock creek; latter langham island. plus pcr negative.
created new script for running test data: dada2_kankakee_test.R
installed cutadapt with `sudo apt install cutadapt`

## 05/22/25
goal: get cutadapt running *

what i did:
i lifted the primers from the paper
here is the command i used to trim the primers off of the first subset fq
```
cutadapt \
-g ACTGGGATTAGATACCCC \
-G TAGAACAGGCTCCTCTAG \
-e 0.25 \
-o Kankakee-1-rep1_S1_L001_R1_001-subset-trimmed.fastq \
-p Kankakee-1-rep1_S1_L001_R2_001-subset-trimmed.fastq \
Kankakee-1-rep1_S1_L001_R1_001-subset.fastq Kankakee-1-rep1_S1_L001_R2_001-subset.fastq
```
both primers were trimmed >9900 times, i.e., most of the reads. that looks good.
i need to add the maximum error rate of 0.25; done, and nearly all reads were trimmed

put all of that in the wrapper. don't know if it works yet. i'm going to test it tomorrow morning

## 05/23/25
it is tomorrow morning

to do:
test cutadapt code block, fix if needed *
add argparse to dada2 code + add function call to wrapper *
figure out best way to manage scripts and output in different directories + hardcoding

what i did:
cutadapt works with a subprocess call and trimms all reads
output writes to a log file
made a dir `src` to store source code
pipeline creates output directory and stores everything there
dada2 runs with trimmed data
added argparse to denoising script
    takes path to trimmed reads, pipeline outdir, and name of asv table to be saved (for later use)
downloaded database from https://doi.org/10.5061/dryad.ghx3ffbz3
cleaned up file structure a bit
    all of dr. picq's data is in subdirectories of test/
    scripts called by the pipeline are in src
    any other functional scripts are in tools

installing dada2 on wsl:
biocmanager is not working
ran `sudo apt install libxt6 libxt-dev` to create the required directory
then in R i ran `install.packages("BiocManager")` and it downloaded
i then ran `BiocManager::install("dada2")`
installed argparse in R

CHANGES TO DADA2:
per the paper, i set truncq to 0 and truncLen to (95, 95)
write asv table out to tsv for next steps?

NOTE:
currently, the wrapper runs entirely in the parent directory for the repo
i have an absolute path for the output dir
i will grab absolute paths to scripts
call scripts from their absolute path + pass input data as an absolute path + grab absolute path for output data to pass to next script

TO-DO:
add argparse to dada2 script and function call to wrapper (pass outdir and trimmed read data) *
this is slow as balls, i need to make it quicker (or maybe i don't honestly)
    problem is loading dada2 and argparse. i don't think i can do anything about the runtime of the tool itself tho
i also need to figure out what info i actually want to write to the log

for tomorrow:
i'm going to try to have the dada2 script write the asv table out fasta format for vsearch
the files are in bmogi/dada2_test/
then i'm going to run this and the database through vsearch

## 05/24/25
to-do:
write asv table to fasta format *
run this file and db through vsearch *
figure out how to map the asvs back to samples? *
    i'm losing some info by converting to fasta
    maybe thats what the metadata is for?
    though i think this is for the stats part later...?

what i did:
denoise_reads.R now writes seqtab.nochim asv table to fasta format using biostrings
edited format_datapase.py to write out taxa map from 1117 db; vsearch needs
    fasta with id + seq
    tsv with id + taxonomy
vserach runs with asv fasta and db fasta. then, hits tsv is merged with db tsv to map asvs to taxa
moved vsearch code to wrapper