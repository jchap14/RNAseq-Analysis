###infer_experiment.py will tell if reads are SE, PE, stranded v unstranded
### usage infer_experiment.py -r REFGENE_BED12 -i $name.Aligned.sortedByCoord.out.bam
### just do on qlogin

module add rseqc
infer_experiment.py \
-r /srv/gsfs0/projects/snyder/chappell/Annotations/GENCODE-v19-GRCh37-hg19/gencode.v19.annotation.bed \
-i ./BAMS_trimmed/O1.Aligned.sortedByCoord.out.bam
