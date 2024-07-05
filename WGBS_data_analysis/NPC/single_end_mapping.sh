#!/usr/bin/bash
        set -u

        keywd=*.u.*
        i=0
        name=()

for line in `ls ${keywd}`
do
        name[${i}]=${line%%"."*}
        let i=${i}+1
done

for sample in   NPC_B_1

do

        echo " processing ${sample}"
        bismark  --parallel 4 --bam --bowtie2 -un  /home/yli/ref/genome/male.mm10/  ${sample}\.fq.gz -o ${sample}\_se
        wait
        cd ${sample}\_se
        deduplicate_bismark -s --bam ${sample}*.bam
        wait

   # bismark_methylation_extractor --bedGraph --gzip --cutoff 2 ${sample}*pe.deduplicated.bam


done
