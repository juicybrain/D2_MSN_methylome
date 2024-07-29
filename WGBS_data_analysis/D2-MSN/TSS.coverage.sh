# computeMatrix scale-regions -S  wgbs_D2_Rep1.bam.s.bw wgbs_D2_Rep2.bam.s.bw  EM_D2_Rep1.bam.s.bw  EM_D2_Rep2.bam.s.bw  -R mm10.90.chr.genes.bed  --numberOfProcessors 4 --regionBodyLength 10000 -a 5000 -b 5000 --binSize 100  -o WGBS_EMseq_genes.scale.5k.coverage.test.gz
	computeMatrix reference-point -S wgbs_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw wgbs_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw -R mm10.90.chr.genes.tss.bed   --binSize 10 -b 2000 -a 2000 -o WGBS_EMseq_tss.test.rp.gz


	plotProfile -m WGBS_EMseq_tss.test.rp.gz  --colors crimson red navy blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup   --refPointLabel TSS  --plotHeight 9  --plotWidth 9 --yMin 0 --yMax 1  --yAxisLabel "meCpG level" -out WGBS_EMseq_TSS_me_CpH.svg --plotTitle "TSS"


# plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --kmeans 4
#  plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar'  --kmeans 4
