#!/bin/bash
#for x in `/bin/ls *.Aligned.sortedByCoord.out.bam` ; do bash ./remove_rRNA_from_TRANSCRIPTOME_bam_2.sh $x; done
jobname=`echo $1`
name=`basename $1 .Aligned.toTranscriptome.out.bam`
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $name.remove_rRNA
#$ -j y
cd .
echo "$name"
echo "Total_Alignments"
samtools view $1 | wc -l
echo "rRNA_Alignments"
samtools view $1 | grep -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` | wc -l
echo "non-rRNA_Alignments"
samtools view $1 | grep -v -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` | wc -l
echo "%_rRNA"
EOF
qsub -V -cwd -l h_vmem=4G -pe shm 8 /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh

##################	Output a report of rRNA % ## uncomment below to generate
# cat remove_rRNA_from_TRANSCRIPTOME_bam.sh.o* | xargs -n8 | tr ' ' '\t' > rRNA_in_transcriptome_alignment.report.txt
# rm remove_rRNA_from_TRANSCRIPTOME_bam.sh.e* remove_rRNA_from_TRANSCRIPTOME_bam.sh.o*


##################	to generate rRNA-free bam for RSEM (*.toTranscriptome.no_rRNA.bam)
# samtools view -h -o $name.sam $1
# cat $name.sam | grep -v -E `cat /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.rRNA-mtRNA_names.txt` > $name.no_rRNA.sam
# samtools view -bS $name.no_rRNA.sam > $name.toTranscriptome.no_rRNA.bam
# rm $name.no_rRNA.sam $name.sam $name.no_rRNA.sam
