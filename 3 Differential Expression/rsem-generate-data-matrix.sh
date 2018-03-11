##### Set the experiment name & run these scripts to get expected counts & TPM matrices
### Run on login node, it works fast

## set experiment name
ExpName=Ajay_HighShear

## make a geneCounts directory & move files there
mkdir geneCountMatrices
mv `find . -name "*.genes.results"` ./geneCountMatrices/
mv `find . -name "*.isoforms.results"` ./geneCountMatrices/
cd geneCountMatrices/

## get expected counts for *.genes.results
module add rsem
rsem-generate-data-matrix `find . -name "*.genes.results" | sort | tr '\n' ' '` > Cnts.tmp
sed 's/.genes.results//g' Cnts.tmp > $ExpName.geneCounts.matrix

## get expected counts for *.isoforms.results
# rsem-generate-data-matrix `find . -name "*.isoforms.results" | sort | tr '\n' ' '` > iCnts.tmp
# sed 's/.isoforms.results//g' iCnts.tmp > $ExpName.isoformCounts.matrix

## get TPMs for *.genes.results
# perl ~/scripts/generate_TPM_matrix.pl `find . -name "*.genes.results" | sort | tr '\n' ' '` > TPMs.tmp
# sed 's/.genes.results//g' TPMs.tmp > $ExpName.geneTPMs.matrix

## get TPMs for *.isoforms.results
# perl ~/scripts/generate_TPM_matrix.pl `find . -name "*.isoforms.results" | sort  | tr '\n' ' '` > iTPMs.tmp
# sed 's/.isoforms.results//g' iTPMs.tmp > $ExpName.isoformTPMs.matrix

## remove tmp files
rm TPMs.tmp Cnts.tmp iCnts.tmp iTPMs.tmp

# Use only if individual results needed, or aren't in directory tree
# ‘rsem-generate-data-matrix’ extracts an expected counts matrix from RSEM expression results
# 
# Usage:
# rsem-generate-data-matrix sampleA.genes.results \ sampleB.genes.results \ ... > output_name.counts.matrix
# 
# module add rsem
# rsem-generate-data-matrix \
# ./$SAMPLE1.genes.results \
# ./$SAMPLE2.genes.results \
# ./$SAMPLEn.genes.results \
# > $NAME.geneCnts.matrix
# 
# perl ~/scripts/generate_TPM_matrix.pl \
# ./RNAseq_SLO_2016.02.22/TRIMMED-ANALYSIS/gene_results/L1.genes.results \
# ./RNAseq_SLO_2016.02.22/TRIMMED-ANALYSIS/gene_results/S1.genes.results \
# ./RNAseq_SLO_2016.02.22/TRIMMED-ANALYSIS/gene_results/O1.genes.results \
# ./RNAseq_SLO_2016.08.18/HiSeq/post_Trim_analysis/S.trim/S.trim.genes.results \
# ./RNAseq_SLO_2016.08.18/HiSeq/post_Trim_analysis/L.trim/L.trim.genes.results \
# ./RNAseq_SLO_2016.08.18/HiSeq/post_Trim_analysis/O.trim/O.trim.genes.results \
# > SLO.2x1701.geneTPMs.matrix
# 
