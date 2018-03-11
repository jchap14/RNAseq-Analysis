#!/bin/bash
# Used to index the genome/annotations for RSEM
# qsub -V -cwd -l h_vmem=4G -pe shm 12 ./RSEM-prepare-hg19-index.sh

rsem-prepare-reference --gtf /srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/Genes/genes.gtf \
/srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/WholeGenomeFasta/genome.fa hg19
