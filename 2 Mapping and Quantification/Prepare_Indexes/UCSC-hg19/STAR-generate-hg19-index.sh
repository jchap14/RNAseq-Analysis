#!/bin/bash

# Preparing genome index for STAR -> this is only required once per reference genome
# qsub -V -cwd -l h_vmem=4G -pe shm 12 ./STAR-generate-hg19-index.sh

cd .
STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir /srv/gsfs0/projects/snyder/chappell/Annotations/STAR_genome_hg19_directory \
--genomeFastaFiles /srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/WholeGenomeFasta/genome.fa \
--sjdbGTFfile /srv/gsfs0/projects/snyder/chappell/Annotations/UCSC-hg19/Genes/genes.gtf \
--sjdbOverhang 100 --outFileNamePrefix STAR_genome_hg19_directory
cd .
