#! /bin/bash

## we use a for loop to run de novo assembly using unicycler tool for each sample

for sample in `cat list.txt`
do
unicycler -1 eanbit_amr/${sample}_trim1.fastq -2 eanbit_amr/${sample}_trim2.fastq -o assembly/${sample}_out

## check the quality of the assembly using quast tool for each sample
quast.py assembly/${sample}_out/assembly.fasta -o quality_assembly
done

