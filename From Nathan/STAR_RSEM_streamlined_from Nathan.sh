mkdir $1;
STAR --genomeDir /mnt/data/annotations/by_release/heterokaryon_GRCh38_GRCm38/STAR_2.4.1b.GRCh38_GRCm38_PhiX_ERCCv2/ \
     --outFileNamePrefix $1/Het_RNASeq_$1_ \
     --readFilesIn /mnt/lab_data/kundaje/projects/heterokaryon/RNASeq/fastq/merged/$1.R1.fastq.gz \
                   /mnt/lab_data/kundaje/projects/heterokaryon/RNASeq/fastq/merged/$1.R2.fastq.gz \
     --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD \
     --outFilterMultimapNmax 20 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --sjdbScore 1 \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM \
     --outWigStrand Stranded \
     --twopassMode Basic \
     --twopass1readsN -1 \
     --runThreadN 16;
  
cd $1;
mv Het_RNASeq_$1_Aligned.toTranscriptome.out.bam Tr.bam; 
cat <( samtools view -H Tr.bam ) <( samtools view -@ 16 Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S 1000000000 -T ./ | tr ' ' '\n' ) | samtools view -@ 16 -bS - > Het_RNASeq_$1_Aligned.toTranscriptome.out.bam;
  
/users/nboley/src/rsem-1.2.22/rsem-calculate-expression --bam --estimate-rspd --no-bam-output --seed 12345 -p 16 --paired-end --forward-prob 0 Het_RNASeq_$1_Aligned.toTranscriptome.out.bam /mnt/data/annotations/indexes/RSEM_indexes/rsem_1.2.22.GRCh38_GRCm38_PhiX_ERCCv2/ quant