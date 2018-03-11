#!/bin/bash

## run command:
## for x in `/bin/ls *.sortedByCoord.out.bam` ; do bash split_rRNA_from_BAM.sh $x; done

### RseQC::split_bam.py will split the original BAM file into 3 small BAM files:
# 1. *.in.bam: reads mapped to exon regions of the gene list
# 2. *.ex.bam: reads that cannot be mapped to regions of the original gene list.
# 3. *.junk.bam: qc-failed or unmapped reads. 
# ^ This removes rRNA reads & unmapped reads from BAM files

##### NOTE: This script can only quantify/remove rRNA alignments from GENOME-aligned BAM (*sortedByCoord.bam)
## genome-aligned BAMs will not work with RSEM-calculate-expression...
## Use a different script to quantify/remove rRNA from transcriptome aligned BAM

## NOTE: unnecessary to remove mitochondrial reads because they are not quantified by RSEM

##### load required modules
module add rseqc

##### set variables
bamfile=`echo $1`
name=`basename $1 .Aligned.sortedByCoord.out.bam`
annoDir="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19"

##### create a tempscript for queue sub
cat > /tmp/tempscript.sh << EOF
#!/bin/bash
#$ -N $name.rRNAsplit
#$ -j y
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -pe shm 12
echo "STARTING split_bam.py"
split_bam.py -i $bamfile -r $annoDir/hg19_rRNA.bed -o $name
EOF
qsub /tmp/tempscript.sh
sleep 1
rm /tmp/tempscript.sh

##################	CLEANUP
# rm *.ex.bam *.in.bam *.junk.bam

##################	Output a report of rRNA % ## uncomment below to generate
# find . -name "rRNA_MT_removal_GENCODE.v19.sh.o*" | xargs -n1 cat > a.txt
# cat a.txt | tr -d ' ' > b.txt
# cat b.txt | xargs -n5 | tr ' ' '\t' > c.txt
# cat c.txt | tr ':' '\t' > rRNA_QC_report.txt
# rm a.txt b.txt c.txt