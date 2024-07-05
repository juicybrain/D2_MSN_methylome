#!/usr/bin/bash

for sample in wgbs_D2_1_S19  wgbs_D2_2_S20 EM_D2_1_S17 EM_D2_2_S18

do

    bismark --bowtie2 --gzip -X 700 -u 31500000 -un --parallel 6 --bam  /home/yli/ref/genome/male.mm10  -1  $sample\_R1_val_1.fq.gz -2 $sample\_R2_val_2.fq.gz  -o ./${sample}\_bismark_bowtie_mm10

    wait
    cd ${sample}\_bismark_bowtie_mm10
                                                                                                 
    wait
    deduplicate_bismark --bam ${sample}*pe.bam
    wait
    samtools sort -@ 10 ${sample}*pe.deduplicated.bam -o ${sample}.pe.deduplicated.s.bam
    picard CollectWgsMetrics I=${sample}.pe.deduplicated.s.bam O=collect_wgs_metrics.txt R=/home/yli/ref/genome/male.mm10/male.mm10.fasta &
    qualimap bamqc -bam ${sample}.pe.deduplicated.s.bam -nw 400 -hm 3 --java-mem-size=64G &
    wait
    bismark_methylation_extractor -p --bedGraph --gzip --comprehensive ${sample}*pe.deduplicated.bam
        cd ..

done
