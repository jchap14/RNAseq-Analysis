#!/bin/bash
# find . -name "*.noMito.bam" | xargs -n1 qsub -V -cwd -l h_vmem=4G -pe shm 12 ./bam2fastq.sh

name=`basename $1 .noMito.bam`

#usage: samtools bam2fq input.bam > output.fastq
echo "samtools bam2fq $1 > $name.no_rRNA.fastq"
samtools bam2fq $1 > $name.no_rRNA.fastq

echo "paste - - - - - - - - < $name.no_rRNA.fastq | tee >(cut -f 1-4 | tr '\t' '\n' > $name.no_rRNA.fastq_1) | cut -f 5-8 | tr '\t' '\n' > $name.no_rRNA.fastq_2"
paste - - - - - - - - < $name.no_rRNA.fastq | tee >(cut -f 1-4 | tr '\t' '\n' > $name.no_rRNA.fastq_1) | cut -f 5-8 | tr '\t' '\n' > $name.no_rRNA.fastq_2

rm $name.no_rRNA.fastq
gzip $name.no_rRNA.fastq_1
gzip $name.no_rRNA.fastq_2