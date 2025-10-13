library(tidyverse)
library(qiime2R)
library(phyloseq)

setwd("C:/Users/bmogi/OneDrive/Documents/UniDocs/MSThesis/physeq")

# importing code adapted from yanxianl on qiime2 forums (https://forum.qiime2.org/t/counting-unique-species-in-r/11771/5)

# Make a phyloseq object in R
metadata <- read_tsv("subset_kankakee_metadata.tsv", comment = "#q2:type") 

table <- read_qza("feature-table.qza")
count_tab <- table$data %>% as.data.frame() 

taxonomy <- read_qza("final_taxa.qza")
tax_tab <- taxonomy$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Consensus)

tree <- read_qza("rooted_tree.qza")

ps <- phyloseq(otu_table(as.matrix(count_tab), taxa_are_rows = T),
               phy_tree(tree$data), 
               tax_table(as.matrix(tax_tab)), 
               sample_data(metadata %>% column_to_rownames("#SampleID")))

ps_no_pcr <- subset_samples(ps, categories != "NA") # drop PCR negative so i can get a better look at the other categories when graphing

asv_richness <- estimate_richness(ps_no_pcr, measures = "Observed") # number of observed ASVs per sample

shannon_div <- estimate_richness(ps_no_pcr, measures = "Shannon") # Shannon diversity per sample

plot_richness(ps_no_pcr, measures=c("Observed", "Shannon"), color="categories") # i can probably play around with this. i'm not really sure what i'm looking for honestly.

shannon_div <- estimate_richness(ps_no_pcr, measures = "Shannon")