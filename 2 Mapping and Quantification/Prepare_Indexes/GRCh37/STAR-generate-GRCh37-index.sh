#!/bin/bash

# Preparing genome index for STAR -> this is only required once per reference genome
# qsub -V -cwd -l h_vmem=4G -pe shm 12 ./STAR-generate-GRCh37-index.sh

cd .
STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir /srv/gsfs0/projects/snyder/chappell/Annotations/STAR_genome_GRCh37_directory \
--genomeFastaFiles /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE/GRCh37.p13.genome.fa \
--sjdbGTFfile /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE/gencode.v19.annotation.gtf \
--sjdbOverhang 100 --outFileNamePrefix STAR_genome_GRCh37_directory
cd .
