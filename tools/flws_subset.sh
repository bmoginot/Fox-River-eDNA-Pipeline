cd data/reads

for i in *
do
    zless $i | head -n 40000 > ../subset/$i
done

cd ../subset

for f in *.fastq.gz
do
    mv "$f" "${f%.fastq.gz}-subset.fq"
done