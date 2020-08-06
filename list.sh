#! /bin/bash

## We used a for loop to create a list of sample names and edited it to remain with only uniq sample IDs. Then remove all the intermediate files

for i in `ls *.fastq`
do
echo $i >> list.txt
done

sed 's\_trim1.fastq\\g' list.txt > list_1.txt ; sed 's\_trim2.fastq\\g' list_1.txt > list_2.txt | sort ; uniq list_2.txt > list_3.txt
rm list.txt ; rm list_1.txt ; rm list_2.txt ; mv list_3.txt list.txt
