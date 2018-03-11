#!/bin/bash

##### get read distribution numbers for BAM files
## for x in `/bin/ls *.sortedByCoord.out.bam` ; do bash get_read_distribution.sh $x; done

##### load required modules
module add samtools
module add rseqc

##### set variables
bamfile=`echo $1`
name=`basename $1 .Aligned.sortedByCoord.out.bam`
annoFile="/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/Genes/genes.bed12"

##### create a tempscript for queue sub
cat > $name.tempscript.sh << EOF
#!/bin/bash
#$ -N $name.distribution
#$ -j y
#$ -V
#$ -cwd
#$ -l h_vmem=8G
#$ -pe shm 6

## get the total number of alignments
read_distribution.py -i $bamfile -r $annoFile > $name.distribution.txt
## remove extra files
# rm
EOF
## submit tempscript & cleanup
qsub $name.tempscript.sh
sleep 1
rm $name.tempscript.sh