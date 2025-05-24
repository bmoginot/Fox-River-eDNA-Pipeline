library(dada2)
library(argparse)
library(Biostrings)

parser <- ArgumentParser(description="denoise trimmed reads and output ASV table")
parser$add_argument("-i", "--input", help="path to trimmed reads", required=TRUE)
parser$add_argument("-o", "--output", help="path to pipeline outdir", required=TRUE)
parser$add_argument("-f", "--fasta", help="what to name the ASV fasta (for vsearch)", required=TRUE)

args <- parser$parse_args()
indir <- args$input
outdir <- args$output
outfasta <- args$fasta

path <- file.path(indir) # locate data
list.files(path)

# divide forward and reverse reads and get sample names
fnFs <- sort(list.files(path, pattern="_R1_001-subset-trimmed.fastq", full.names = TRUE)) # check if my file names are in the same format
fnRs <- sort(list.files(path, pattern="_R2_001-subset-trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) # uses firsts field of fq file name

# get quality of reads
# plotQualityProfile(fnFs[1:2])
# plotQualityProfile(fnRs[1:2])
# quality falls off precipitously around the same position for the forward and reverse reads. reverse reads however are worse earlier, expectedly.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(outdir, "filtered_reads", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, "filtered_reads", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(95,95),
                     maxN=0, maxEE=c(2,2), truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)
# lost a few hundred reads, everything seems to be in order here

# learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
# yeah i think these look fine

# sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]] # dada class objects contain many diagnostics that can be interrogated further

# merge
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]]) # merge objects are lists of dataframes with information about each merge. make sure most reads successfully merge

# construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # how many asvs are there?
table(nchar(getSequences(seqtab)))
# there are a few sequences that are really long but i think they were removed as chimeras in the next step

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# still have most of our reads

# look at number of reads that made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
# looks good to me; only about 1000 reads lost for each sample

# some code to convert the asv table to fasta format for vsearch
asv_seqs <- colnames(seqtab.nochim)
asv_fasta <- DNAStringSet(asv_seqs)
names(asv_fasta) <- paste0("ASV_", seq_along(asv_seqs))
write_asv_out <- file.path(outdir, outfasta)
writeXStringSet(asv_fasta, write_asv_out)