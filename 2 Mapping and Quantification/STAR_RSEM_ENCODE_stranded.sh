#!/bin/bash

##### STAR mapping / RSEM quantification pipeline #####

## run command:
## for x in `/bin/ls *.trim.R1.fq.gz` ; do bash STAR_RSEM_ENCODE_stranded.sh $x; done

## set variable names
read1=`echo $1` #gzipped fastq file for read1
NAME=`basename $1 .trim.R1.fq.gz` #for trimmed
read2=$NAME.trim.R2.fq.gz #gzipped fq file for read2, use "" if single-end

## add required modules
module add STAR
module add ucsc_tools
# module add rsem     #in conda env
# module add samtools #in conda env
# module add r        #in conda env

# STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
STARgenomeDir="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/STAR_genome_GRCh37_directory"
RSEMrefDir="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/RSEM_genome_GRCh37_directory/GRCh37"
nThreadsSTAR="12" # number of threads for STAR
nThreadsRSEM="12" # number of threads for RSEM

## put the STAR & RSEM commands here ##
cat > $NAME.tempscript.sh << EOF
#!/bin/bash -l
#SBATCH --job-name $NAME.STAR_RSEM
#SBATCH --output=$NAME.STAR_RSEM.out
#SBATCH --mail-user jchap14@stanford.edu
#SBATCH --mail-type=ALL
# Request run time & memory
#SBATCH --time=5:59:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --export=ALL
#SBATCH --account=mpsnyder

mkdir $NAME
cd $NAME

###### STAR command
echo "STARTING STAR"
STAR --genomeDir $STARgenomeDir  --readFilesIn ../$read1 ../$read2 \
 --outFileNamePrefix $NAME. --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20 \
 --alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat \
 --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000 \
 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM  $STARparStrand

###### bedGraph generation, now decoupled from STAR alignment step
echo "STARTING STAR bedGraph generation"
mkdir Signal
STAR --runMode inputAlignmentsFromBAM --inputBAMfile $NAME.Aligned.sortedByCoord.out.bam \
--outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix ./Signal/ \
--outWigReferencesPrefix chr

# move the signal files from the subdirectory
echo "move the signal files from the subdirectory"
mv Signal/Signal*bg .
# 
# 
# ###### bigWig conversion commands
echo "Start bigWig conversion"
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

# stranded data
str[1]=-; str[2]=+;
for istr in 1 2
do
for imult in Unique UniqueMultiple
do
	grep ^chr Signal.\$imult.str\$istr.out.bg > sig.tmp
	bedSort sig.tmp sig.tmp
bedGraphToBigWig sig.tmp  chrNL.txt Signal.\$imult.strand\${str[istr]}.bw
done
done
rm sig.tmp
##########################################################################################
######### RSEM
echo "Prepare for RSEM: sorting BAMs"

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
mv $NAME.Aligned.toTranscriptome.out.bam Tr.bam 


# paired-end data, merge mates into one line before sorting, and un-merge after sorting
# paired-end data, merge mates into one line before sorting, and un-merge after sorting
echo "cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf \"%s\", \$0 \" \"; getline; print}' | sort -S 60G -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > $NAME.Aligned.toTranscriptome.out.bam"
cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf "%s", \$0 " "; getline; print}' | sort -S 60G -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > $NAME.Aligned.toTranscriptome.out.bam
'rm' Tr.bam

###### RSEM command
echo "STARTING RSEM"
rsem-calculate-expression --bam --estimate-rspd --no-bam-output --seed 12345 -p $nThreadsRSEM --paired-end \
--forward-prob 0 $NAME.Aligned.toTranscriptome.out.bam $RSEMrefDir $NAME >& $NAME.Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file $NAME.pdf, which contains multiple plots
echo "STARTING RSEM-plot-model"
echo rsem-plot-model $NAME $NAME.pdf
rsem-plot-model $NAME $NAME.pdf

EOF
# qsub then remove the tempscript
sbatch $NAME.tempscript.sh 
sleep 1
rm $NAME.tempscript.sh
