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
module add samtools
module add rsem
module add r
module add ucsc_tools

# STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
STARgenomeDir="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/STAR_genome_GRCh37_directory"
RSEMrefDir="/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/RSEM_genome_GRCh37_directory/GRCh37"
nThreadsSTAR="12" # number of threads for STAR
nThreadsRSEM="12" # number of threads for RSEM

# executables
STAR=STAR                             
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig              

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn ../$read1 ../$read2 \
 --outFileNamePrefix $NAME --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20 \
 --alignIntronMax 1000000   --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"

# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 
#OPTION: stranded data
STARparStrand=""
STARparWig="--outWigStrand Stranded"

# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM "

# RSEM parameters:
#OPTION: stranded paired end
RSEMparType="--paired-end --forward-prob 0"
      
## put the STAR & RSEM commands here ##
cat > $NAME.tempscript.sh << EOF
#!/bin/bash -l
#SBATCH --job-name $NAME.STAR_RSEM
#SBATCH --output=$NAME.STAR_RSEM.out
#SBATCH --mail-user jchap14@stanford.edu
#SBATCH --mail-type=ALL
# Request run time & memory
#SBATCH --time=5:59:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --export=ALL
#SBATCH --account=mpsnyder

mkdir $NAME
cd $NAME

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

###### STAR command
echo "STARTING STAR"
echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand
$STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand

###### bedGraph generation, now decoupled from STAR alignment step
# working subdirectory for this STAR run
echo "STARTING STAR bedGraph generation"
mkdir Signal

echo $STAR --runMode inputAlignmentsFromBAM --inputBAMfile $NAME.Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
$STAR --runMode inputAlignmentsFromBAM --inputBAMfile $NAME.Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

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
$bedGraphToBigWig sig.tmp  chrNL.txt Signal.\$imult.strand\${str[istr]}.bw
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
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $NAME.Aligned.toTranscriptome.out.bam $RSEMrefDir $NAME >& $NAME.Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType $NAME.Aligned.toTranscriptome.out.bam $RSEMrefDir $NAME >& $NAME.Log.rsem

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
