#!/usr/bin/bash

wd="/home/yli/fenglab/YuxiangTempData/41_WGBS"

for sample in wgbs_D2_1_S19  wgbs_D2_2_S20 EM_D2_1_S17 EM_D2_2_S18

do

    bismark   --bowtie2 --gzip --parallel 4  --bam  /home/yli/ref/genome/lambdaDNA -1 $sample\_R1_val_1.fq.gz -2 $sample\_R2_val_2.fq.gz -o ./${sample}_bismark_bowtie_lambda

    wait
    cd ${sample}_bismark_bowtie_lambda
    deduplicate_bismark --bam ${sample}*.bam
    wait
    bismark_methylation_extractor -p --bedGraph --gzip --CX --ignore_3prime 3 --ignore_3prime_r2 3 --cytosine_report --genome_folder ~/ref/genome/lambdaDNA/  *.deduplicated.bam

    zcat *.CX_report.txt.gz | grep $'\t'CG[ATCG]$ > CG.lambda.cov
    zcat *.CX_report.txt.gz | grep $'\t'CC[ATCG]$ > CC.lambda.cov
    zcat *.CX_report.txt.gz | grep $'\t'CT[ATCG]$ > CT.lambda.cov
    zcat *.CX_report.txt.gz | grep $'\t'CA[ATCG]$ > CA.lambda.cov

    cat CA.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CA", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
    cat CC.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CC", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
    cat CT.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CT", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
    cat CG.lambda.cov  |awk 'BEGIN{me=0;unme=0;OFS="\t"}{me=me+$4;unme=unme+$5}END{print "CG", me, unme,1-me/(me+unme)}' >> coversion_rate.CX.txt
    cd ..

done
