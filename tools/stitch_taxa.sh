# very simple command to create aggregate taxa file for phyloseq
cp output/vsearch/retained_vsearch_taxa.tsv taxa.tsv # duplicate vsearch output
tail -n +2 output/bayes/retained_bayes_taxa.tsv >> taxa.tsv # concatenate all but header from bayes output