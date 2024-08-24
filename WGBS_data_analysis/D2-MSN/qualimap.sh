for spl in    wgbs_D2_Rep1.R1_bismark_bt2_pe.deduplicated wgbs_D2_Rep1.R2_bismark_bt2_pe.deduplicated EM_D2_Rep1.R1_bismark_bt2_pe.deduplicated EM_D2_Rep2.R1_bismark_bt2_pe.deduplicated


do

# qualimap bamqc -bam ${spl}\.s.bam -nw 400 -hm 3 --java-mem-size=32G
# picard CollectWgsMetrics I=${spl}.s.bam O=${spl}\collect_wgs_metrics.txt R=/home/yli/ref/genome/male.mm10/male.mm10.fasta
picard CollectGcBiasMetrics I=${spl}\.s.bam  O=${spl}\.gc_bias.txt CHART=${spl}\.gc.pdf S=${spl}\.summary.txt R=/home/yli/ref/genome/male.mm10/male.mm10.fasta
 done

