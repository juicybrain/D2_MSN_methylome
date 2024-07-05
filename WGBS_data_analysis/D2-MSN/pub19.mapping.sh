#!/usr/bin/bash

for sample in   wgbs_D2_1_S19  wgbs_D2_2_S20 EM_D2_1_S17 EM_D2_2_S18

do
    	bismark   --bowtie2  --score_min L,0,−0.2 --gzip --parallel 4 --bam  /home/yli/ref/genome/puc19 -1  $sample\_R1_val_1.fq.gz -2 $sample\_R2_val_2.fq.gz  -o ./${sample}_bismark_bowtie_puc19
        wait
        cd ${sample}_bismark_bowtie_puc19
        deduplicate_bismark --bam ${sample}*.bam
        wait
         bismark_methylation_extractor -p --bedGraph --cutoff 1 --ignore_3prime 3 --ignore_3prime_r2 3 ${sample}*pe.deduplicated.bam
        cd ..

done
