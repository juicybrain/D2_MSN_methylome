#!/usr/bin/bash
wd='/home/fenglab/Projects/26_nucleicRNAseq_OxBS_practice/Fastq2_Feng052620/'

for sample in wgbs_D2_1_S19  wgbs_D2_2_S20 EM_D2_1_S17 EM_D2_2_S18

do

    bismark   --bowtie2 --gzip --parallel 4 --bam  /home/yli/ref/genome/5hmc_std_bowtie_index/  -1  $sample\_R1_val_1.fq.gz -2 $sample\_R2_val_2.fq.gz -o ./${sample}_bismark_bowtie_5hmc_bam
    wait
    cd ${sample}_bismark_bowtie_5hmc_bam
    deduplicate_bismark --bam ${sample}*.bam
    wait
    bismark_methylation_extractor -p --bedGraph --gzip  --ignore_3prime 3 --ignore_3prime_r2 3 ${sample}*pe.deduplicated.bam

    cd ..

done
