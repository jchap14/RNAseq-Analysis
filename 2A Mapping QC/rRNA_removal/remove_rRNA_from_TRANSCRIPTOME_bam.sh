#!/bin/bash

### NOTE: This can quantify/remove rRNA alignments from TRANSCRIPTOME-aligned BAM (*toTranscriptome.bam)
### This should work with RSEM-calculate-expression...
### Use other script to quantify/remove rRNA from GENOME-aligned BAM (*sorted.byCoordinate.bam)

#RUN COMMAND below ->
# find . -name "*.toTranscriptome.out.bam" | xargs -n1 qsub -V -cwd -l h_vmem=4G -pe shm 4 ./remove_rRNA_from_TRANSCRIPTOME_bam.sh

name=`basename $1 .Aligned.toTranscriptome.out.bam`

##################	Outputs rRNA stats to *.o files  
echo "$name"
echo "Total_Alignments"
samtools view $1 | wc -l
echo "rRNA_Alignments"
samtools view $1 | grep -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` | wc -l
echo "non-rRNA_Alignments"
samtools view $1 | grep -v -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` | wc -l
echo "%_rRNA"

##################	Output a report of rRNA % ## uncomment below to generate
# cat remove_rRNA_from_TRANSCRIPTOME_bam.sh.o* | xargs -n8 | tr ' ' '\t' > rRNA_in_transcriptome_alignment.report.txt
# rm remove_rRNA_from_TRANSCRIPTOME_bam.sh.e* remove_rRNA_from_TRANSCRIPTOME_bam.sh.o*


##################	to generate rRNA-free bam for RSEM (*.toTranscriptome.no_rRNA.bam)
# samtools view -h -o $name.sam $1
# cat $name.sam | grep -v -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` > $name.no_rRNA.sam
# samtools view -bS $name.no_rRNA.sam > $name.toTranscriptome.no_rRNA.bam
# rm $name.no_rRNA.sam $name.sam $name.no_rRNA.sam
