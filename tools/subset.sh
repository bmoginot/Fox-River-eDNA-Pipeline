# subsets kankakee data to 10k reads and renames them to *-subset.fq

for i in KAN*
do
    cd $i
    for j in $(ls)
    do
        zless $j | head -n 40000 > ../$j
        $n 
    done
    cd ..
done

for f in *.fastq.gz
do
    mv "$f" "${f%.fastq.gz}-subset.fastq"
done