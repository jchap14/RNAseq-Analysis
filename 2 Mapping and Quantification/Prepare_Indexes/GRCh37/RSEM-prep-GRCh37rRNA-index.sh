#!/bin/bash
# Used to index the genome/annotations for RSEM
# qsub -V -cwd -l h_vmem=4G -pe shm 12 ./RSEM-prep-GRCh37rRNA-index.sh

rsem-prepare-reference --gtf /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE/gencode.v19.rRNA.annotation.gtf \
/srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE/GRCh37.p13.genome.fa GRCh37rRNA
