#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l p=3
#$ -l vf=40G


path='Your/Working/Directory'
input="1-RPF_R2.fq.gz"    ###### We only use R2 of each ribo-seq data
index="new_sample-1_RPF" ##### output name

date 
echo 'trim-1-doing'
cutadapt -g GTGCTCTTCCGATCT -g TGTGCTCTTCCGATC -g GTGTGCTCTTCCGAT -g CGTGTGCTCTTCCGA \
-g ACGTGTGCTCTTCCG -g GACGTGTGCTCTTCC -a AGATCGGAAGAGCGT -a GATCGGAAGAGCGTC -a ATCGGAAGAGCGTCG \
-a TCGGAAGAGCGTCGT -a CGGAAGAGCGTCGTG -a GGAAGAGCGTCGTGT -e 0.1 -O 5 -n 2 -f fastq -m 11  \
--match-read-wildcards -o $index-trim-step1.fq.gz $path/$input > $index-trim-step1-adaptor.metrics \
--info-file $index-trim-step1.log

## rm 10nt random sequence
echo 'trim-2-doing'
python modified_demux_single_end.py \
--fastq_2 $index-trim-step1.fq.gz --out_file_2 $index-trim-step2.fq.gz --length 10


## rm low quality reads
java -jar trimmomatic-0.36.jar SE \
-phred33 $index-trim-step2.fq.gz $index-trim-step3.fq.gz SLIDINGWINDOW:4:15 MAXLEN:50 MINLEN:25 \
MAXINFO:30:0.9 -trimlog $index-trim-step3.log -threads 3

## rm rRNA
bowtie2 --un $index-unrRNA.fq -q --sensitive-local -p 3 -N 0 -L 20 \
-x ref_data/genome/hg38-rRNA/h38-rRNA \
-U $index-trim-step3.fq.gz -S $index-rRNA.sam 2> $index-rRNA.log

## rm mitochondrial RNA and tRNA
echo 'trim Mitochondrial'
bowtie2 --un $index-unMT.fq -q --sensitive-local -p 3 -N 0 -L 20 \
-x /ref_data/hg38-Mitochondrial/hg38-Mitochondrial \
-U $index-unrRNA.fq -S $index-MT.sam 2> $index-MT.log

## rm tRNA
bowtie2 --un $index-untRNA.fq -q --sensitive-local -p 3 -N 0 -L 20 \
-x ref_data/hg38-tRNA/hg38-tRNA \
-U $index-unMT.fq -S $index-tRNA.sam 2> $index-tRNA.log


echo $sample 'Head 3 base cutting'
java -jar trimmomatic-0.36.jar SE \
-phred33 $index-untRNA.fq $index-untRNA-cuthead3.fq \
HEADCROP:3 -trimlog $index-trim-step3.log -threads 3


### mapping to genome
STAR --runMode alignReads --runThreadN 3 --genomeDir ref_data/genome/genecode/gencode_release36_GRCh38.p13/STAR_index \
--readFilesIn $index-untRNA-cuthead3.fq --outSAMunmapped Within \
--outFilterMultimapNmax 1 --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2 \
--outSAMattributes All  --outSAMtype BAM Unsorted --outFilterType BySJout \
--outReadsUnmapped Fastx --outFileNamePrefix $index


### extract unique reads
samtools view -H *.bam > 0-head.txt
samtools view *.bam | grep 'NH:i:1[^0-9]' > 0-all.txt
cat 0-head.txt 0-all.txt > 0-all.sam
samtools view -bS 0-all.sam > 0-all.bam
samtools sort -@ 3 0-all.bam $index-sorted.bam #这一句的参数要注意
samtools index $index-sorted.bam.bam
rm 0-head.txt
rm 0-all.txt
rm 0-all.sam
rm 0-all.bam


### remove PCR duplicates
python new.py --bam $index-sorted.bam.bam --out_file $index.rmDup.bam \
--metrics_file out.rmDup.metrics

less *.rmDup.metrics|perl -e 'my $total;my $rm;print "sample\ttotal_No\trm_No\tusable_reads\n";while(<>){chomp;if($_=~/^randomer/){}else{my @a=split/\t/;$total+=$a[1];$rm+=$a[2];}}$usable=$total-$rm;print "*\t$total\t$rm\t$usable\n";'>count_usable.txt
samtools sort -@ 3 -m 6G -f $index.rmDup.bam  $index.rmDup.sorted.bam
samtools index $index.rmDup.sorted.bam 
rm $index.rmDup.bam
 
#### Count read number
featureCounts -s 1 -t CDS -g gene_id -T 3 -a gencode.v36.annotation.gtf -o ${index}_counts.txt $index.rmDup.sorted.bam

