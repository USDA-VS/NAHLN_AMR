#!/usr/bin/env bash

fasta=*.fasta
name=`echo $fasta | sed 's/.fasta//'`
grep ">" $fasta >> contig.txt
length=0
for i in `cat contig.txt`; do
    num=`echo $i | sed -r 's/(.*?)_length_//' | sed -r 's/_cov_(.*?)//' `
    if [ $num -gt 1000 ]; then
    length=$((length + num))
    fi
done
rm contig.txt
echo $length >> ${name}_lk.txt
