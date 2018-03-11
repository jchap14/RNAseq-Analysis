#!/bin/bash
# Quantifies transcripts from BAM files 
name=`basename $1 .Aligned.sortedByCoord.out.bam`
# for x in `/bin/ls *.bam` ; do bash ./Cufflinks.sh $x; done
name=`basename $1 .Aligned.sortedByCoord.out.bam`
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $name.cufflinks
#$ -j y
echo "STARTING CUFFLINKS"
cufflinks -p 16 -o $name \
-g /home/jchap14/hi_quota_folder/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.annotation.gtf \
--library-type fr-firststrand $1
EOF
qsub -V -cwd -l h_vmem=2G -pe shm 16 -l h_rt=99:99:99 /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh

##########################################################################################
##########################################################################################
################################ CUFFMERGE -> merge all cufflinks transcripts into one GTF
# find . -name "*_transcripts.gtf" > assemblies.txt



##########################################################################################
##########################################################################################
################################ CUFFDIFF -> merge all cufflinks transcripts into one GTF
# find . -name "*_transcripts.gtf" > assemblies.txt
