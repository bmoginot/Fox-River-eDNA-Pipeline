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
moved all test data from test/ to data/
moved vsearch code to wrapper. seems to be producing expected output.
the paper cites orourke 2021, 2022 for the "hybrid approach" of taxa classification. i'll read through these later to see exactly what they did. may have to rework vsearch, we will see.
added subproccess call for vsearch.
also changed all subprocess calls to also write stderr to the log.
    i need to figure out what i want to write out to log and to console. though, inevitably, all of this info will just be for me.
    maybe tell me when certain things are finished, etc. + write out results to log??

for tomorrow (or later):
read orourke 2021, 2022 and make changes as appropriate
figure out the log
rescript

check these out:
‘qiime rescript evaluate-fit-classifier’
‘qiime feature-classifier classify-sklearn’

## 05/25/25
what i did:
updated vsearch parameters to more closely match qiime2 defaults
qiime2 computes LCA when more than one sequence matches a given ASV. this doesn't seem too hard, i can do this on my own.
    done
at this point, we keep anything classified to family or above and move the rest into the naiive bayesian classifier next.
this is rapidly exiting the scope of my project. i need to download qiime2 and do everything per the paper's methods.

for tomorrow:
scrap everything and use qiime2 instead

## 05/26/25
made a new wrapper to run qiime2
before i can do anything, i need to turn my fastq directory into a qiime archive (.qza)
    in order to do this i need a manifest file
gzip-ed all of the subset reads using `gzip -r data/subset_reads/`
added code to create manifest file, successfully imported reads
cutadapt works as well
    verified by
    unzip output/reads.qza
    `zless | head` for the first sample in data/
    comparing it against the first sample in trimmed_reads.qza after unzipping that as well
dada2 works
for vsearch, i need the asv-sequences-0.qza from dada2. i also need to read in a tsv file and a fasta file for reference, this is done via tools/qiime_format_database.py
vsearch works and outputs tsv and blast6 files
updated readme to include information about files going in and out at each step

## 05/27/25
to-do:
naive bayes classification
    need to train the model with a database via rescript
    then classify sequences using model

what i did:
updated wrapper to print out the current step so i'm not just staring at a blank screen for 5 minutes
trained model with rescript and saved it to rescript_classifer.qza
used model to classify asv sequences and output nb_classification.qza
added multithreading to everything i can (frick dada2)
got rid of log write-out since most output is supressed by qiime anyway (write my own log later)

there's one more thing i wanted to do but i can't remember

## 05/28/25
to-do:
filter vsearch output *
pass only non-family-level seqs to bayes
filter bayes output

what i did:
script (src/parse_vsearch_output.py) that extracts vsearch output, filters for family-level classification, then splits the tsv into seqs that are sufficiently classified and those that need to go to the bayes classifier
integrated this script into naive bayes script

next steps:
filter dada2 asvs using unassigned feature ids from vsearch
import this and send it to the bayesian classifier
do the same for bayes output

## 05/31/25
what i did:
parse_output now processes and filters both vsearch and bayes tsvs
moved all of the functions over to the wrapper

next steps:
i passed the retained vserach tsv to the bayesian classifier, but what i need to do is take the unassigned tsv, get all of the feature ids and use that to grab sequences from dada2 output. i will write a script to do that tomorrow.

## 06/02/25
here's the plan:
take asv fasta output from dada2 and build hashmap for id --> sequence
use feature id column from tsv for seqs that were *unclassified* by vsearch
iterate this vector, index the hashmap, retrieve seqs, rebuild fasta --> bayes classifier
this is so much shit fr fr
i'm going to build this in the wrapper i think that makes more sense

the wrapper now:
takes the qiime archive generated by vsearch.
extracts the taxonomy tsv from it.
splits it by family-level and sub-family-level classification into two tsvs.
makes a hash map from dada2 fasta output.
uses the features ids from the sub-family-level classifications to index that hashmap and rebuild a fasta file of only sequences that were left unclassified by vsearch.
imports this fasta as qza and passes it to the bayesian classifier.
performs the same process on the bayes output.

holy shit everything works so far
everything seems to be in order but i feel like i need a way to check my work. maybe i'll mention this to yoel when he's back.

i also really really need to clean up the output folder and figure out a log

next step: blast with 98% identity and 95% coverage (how to do this in house?)

## 06/04/25
i think i'm just going to do some housekeeping stuff today

what i did:
added some docstrings
created subdirs for dada2, vsearch, and bayes. the parse_output and map_seqs functions also write to vsearch and bayes depending on which input they are called with.
also updated the readme
blasted the unclassified bayes seqs. there were only two that fit the criteria (95% coverage, 98% identity) and they were both human

## 06/10/25
ope its been a minute

to-do:
keep working on documentation
stats...?

what i did:
documentation in the README is up to date (through taxonomic classification)

for tomorrow:
extract biom file from qiime dada2 output and import into phyloseq. import qiime sample data (metadata from dr. picq i believe). merge. go on my merry way

## 06/11/25
i found an R package (qiime2R) that can create a phyloseq object from qiime archives
i need a feature table, a taxonomy table, a tree, and metadata
feature table <-- dada2
metadata <-- in hand from FMNH

need to stitch taxonomy tsvs together (need to do same for fastas)