##### These scripts are for RNAseq QC after BAM_indexing (samtools_index_stats.sh)
## run from ./BAMs/ , but it outputs the FQ_BAM_stats.QC file in ..

##### Get stats for # of (PE) reads in FASTQs & the mapped read %
find .. -name "*.Log.final.out" | sort | sed 's/.Log.final.out//g' > a.txt
find .. -name "*.Log.final.out" | sort | xargs -n1 cat | grep "input reads" | tr -d ' ' | cut -f2 > b.txt
find .. -name "*.Log.final.out" | sort | xargs -n1 cat | grep "Uniquely mapped reads %" | tr -d ' ' | cut -f2 > d.txt
paste -d "\t" a.txt b.txt > c.txt
paste -d "\t" c.txt d.txt | sort > fastq_stats.txt
## add column header to file
echo -e "ID\tFQread#\tUniqMap%" > header
cat header | cat - fastq_stats.txt > temp && mv temp fastq_stats.txt
rm a.txt b.txt c.txt d.txt header

##### Get total alignment number & mitochondrial alignment number
cat bam_stats.txt | xargs -n3 -d "\n" | sort > a.txt
cat a.txt | sed 's/ + 0 in total (QC-passed reads + QC-failed reads)//g' > b.txt
cat b.txt | sed 's/.Aligned.sortedByCoord.out.bam//g' > e.txt
cat e.txt | tr ' ' '\t' | cut -f1,2,5  > Alignment_stats.txt

## add column header to file
echo -e "ID\tBAMalign#\tchrMalign#" > header
cat header | cat - Alignment_stats.txt > temp && mv temp Alignment_stats.txt
rm a.txt b.txt header e.txt

##### Get stats for rRNA read # of (PE) reads in FASTQs
find *.IndexAndStat.o* | sort | cut -f1 -d "." > a.txt
find *.IndexAndStat.o* | sort | xargs -n1 cat | grep ".in.bam" | cut -f2 -d ":" | tr -d ' ' > b.txt
paste -d "\t" a.txt b.txt > rRNA_readNum.txt
## add column header to file
echo -e "ID\trRNAread#" > header
cat header | cat - rRNA_readNum.txt > temp && mv temp rRNA_readNum.txt
rm a.txt b.txt header

############## Combine into 1 report
paste -d "\t" fastq_stats.txt Alignment_stats.txt > a.txt
paste -d "\t" a.txt rRNA_readNum.txt > ../FQ_BAM_stats.QC

## remove extra files
rm fastq_stats.txt Alignment_stats.txt rRNA_readNum.txt a.txt

##### Get duplication stats
## generate Duplication_stats.txt after running Picard-MarkDuplicates.sh
rm Picard-MarkDuplicates.sh.e* Picard-MarkDuplicates.sh.o*
find . -name "*.dedup_metrics.txt" | sort > a.txt
find . -name "*.dedup_metrics.txt" | sort | xargs -n1 cat | grep -e "Unknown" | tr -d ' ' > b.txt
paste -d "\t" a.txt b.txt > c.txt
echo -e "ID\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE" > header 
cat header | cat - c.txt > temp && mv temp Duplication_stats.txt
rm a.txt b.txt c.txt header
