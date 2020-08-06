#! /bin/bash

## Here, we check the quality of each sample using fastqc tool and then get a merged report using multiqc tool

for i in `ls eanbit_amr/*.fastq`
do
fastqc -o quality $i
done

multiqc quality 

