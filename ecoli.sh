
source activate eTrack

path="/home/jkabahita/mproject"
mkdir -p ${path}/results/trim-results
mkdir -p ${path}/results/trim-results/refastqc 
mkdir -p ${path}/results/sam-results
mkdir -p ${path}/results/unicycler-results
mkdir -p ${path}/results/snippy-results
mkdir -p ${path}/results/mlst
mkdir -p ${path}/results/prokka
mkdir -p ${path}/results/abricate
mkdir -p ${path}/results/roary
mkdir -p ${path}/results/iqtree

path2="/home/jkabahita/mproject/results"

## Getting sample ids

cd ~/mproject/fastq

for i in ` ls *trim1*`
do
    basename -s _trim1.fastq ${i} >> ~/mproject/list.txt

    mv ${i}_trim1.fastq.gz ${i}_R1.fastq.gz

    mv ${i}_trim2.fastq.gz ${i}_R2.fastq.gz
done 

cd ~/mproject 

reads="${path}/trim-results"



bwa index -p ecoliINDEX -a is /home/jkabahita/mproject/ecoli-ref.fasta
 
 for i in `cat ${path}/try.txt`
 do
     echo "now running trimmomatic"
     trimmomatic PE  ${path}/fastq/${i}_R1.fastq.gz ${path}/fastq/${i}_R2.fastq.gz \
		  ${path}/results/trim-results/${i}-processed_1P.fastq.gz \
		  ${path}/results/trim-results/${i}-processed_1U.fastq.gz\
		  ${path}/results/trim-results/${i}-processed_2P.fastq.gz\
		  ${path}/results/trim-results/${i}-processed_2U.fastq.gz\
		 -threads 20\
		 -trimlog trim.log \
		 ILLUMINACLIP:/home/jkabahita/Sequencing_adaptors.fasta:2:30:10\
		 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
     #fastqc -o ${path}/results/trim-results/refastqc  ${path2}/trim-results/${i}-processed.fastq.gz


     bwa mem  ecoliINDEX -t 20 \
	${path2}/trim-results/${i}-processed_1P.fastq.gz\
	${path2}/trim-results/${i}-processed_2P.fastq.gz > ${path2}/sam-results/${i}.sam

    samtools view -h -b ${path2}/sam-results/${i}.sam -o ${path2}/sam-results/${i}.bam -@20

    samtools sort ${path2}/sam-results/${i}.bam -o ${path2}/sam-results/${i}.bam -@20

    samtools index ${path2}/sam-results/${i}.bam -@20

    mkdir -p ${path2}/unicycler-results/${i}

    unicycler -1  ${path2}/trim-results/${i}-processed_1P.fastq.gz \
	   -2  ${path2}/trim-results/${i}-processed_2P.fastq.gz \
	   -o  ${path2}/unicycler-results/${i} \
	   --keep 0 \
	   -t 20

   # snippy --cpus 20 --reference  /home/jkabahita/mproject/ecoli-ref.gb\
   #	   --R1 ${path2}/trim-results/${i}-processed_1P.fastq.gz \
#	   --R2 ${path2}/trim-results/${i}-processed_2P.fastq.gz \
#	   --force \
#	   --outdir ${path2}/snippy-results/${i} \
#	   --prefix ${i} \
    #	   --cleanup

    prokka  --force --outdir  ${path}/results/prokka \
	    --prefix ${i} --cpus 20 \
	    ${path2}/unicycler-results/${i}/assembly.fasta

    mlst --threads 20 --scheme ecoli \
	 ${path2}/unicycler-results/${i}/assembly.fasta  | echo -e  "${i} -\t `grep ecoli`" >> ${path2}/mlst/mlst.txt
    
    cut -f1,3-  ${path2}/mlst/mlst.txt > ${path2}/mlst/mlstfinal.txt

    abricate "${path2}/unicycler-results/${i}/assembly.fasta"  > "${path2}/abricate/${i}.txt"

    cat ${path2}/abricate/${i}.txt | cut -f6,10,11,14,15 > ${path2}/abricate/${i}_1.txt

    echo " " >> ${path2}/abricate/abr_final.txt

    echo -e "abricate results for ${i}"  >> ${path2}/abricate/abr_final.txt

    echo " " >> ${path2}/abricate/abr_final.txt

    cat ${path2}/abricate/${i}_1.txt >> ${path2}/abricate/abr_final.txt


 done


  
