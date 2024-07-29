computeMatrix scale-regions -S  wgbs_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw wgbs_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw  -R mm10.90.chr.genes.bed  --numberOfProcessors 4 --regionBodyLength 10000 -a 5000 -b 5000 --binSize 100  -o WGBS_EMseq_genes.scale.5k.me.test.gz



    plotProfile -m WGBS_EMseq_genes.scale.5k.me.test.gz  --colors crimson red navy blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup    --startLabel TSS  --endLabel TES  --plotHeight 9  --plotWidth 12 --yMin 0 --yMax 1  --yAxisLabel "meCpG level" -out WGBS_EMseq_gene_body_me_CpH.svg

#    plotProfile -m WGBS_EMseq_enhancer.scale.5k.me.test.gz  --colors crimson red navy blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup  --startLabel enhancer --endLabel enhancer -out WGBS_EMseq_enhancers_me_CpG.svg --plotHeight 9  --plotWidth 9 --yMin 0 --yMax 0.83  --yAxisLabel "meC level" &

# plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --kmeans 4

#  plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar'  --kmeans 4



computeMatrix scale-regions -S  wgbs_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw wgbs_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep1.CH.cov.gz.GT5.bedGraph.bw EM_D2_Rep2.CH.cov.gz.GT5.bedGraph.bw  -R mm10.90.chr.genes.bed  --numberOfProcessors 4 --regionBodyLength 10000 -a 5000 -b 5000 --binSize 100  -o WGBS_EMseq_genes.scale.5k.me.test.gz

#  plotProfile -m WGBS_EMseq_Str_EnHr.rp.gz --perGroup  --colors red red blue blue  --plotHeight 21  --plotWidth 40 -out WGBS_EMseq_EnHr.svg
#  plotProfile -m WGBS_EMseq_genes.scale.5k.me.test.gz  --colors red red blue blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup  --startLabel TSS --endLabel TES --plotHeight 8  --plotWidth 10 -out WGBS_EMseq_genebody_me.svg

#    plotProfile -m WGBS_EMseq_genes.scale.5k.me.test.gz  --colors crimson red navy blue  --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup  --startLabel TSS --endLabel TES --plotHeight 9  --plotWidth 9 -out WGBS_EMseq_genes_me_CpG.svg  --yMin 0 --yMax 1  --yAxisLabel "meC level" &
    plotProfile -m WGBS_EMseq_genes.scale.5k.me.test.gz  --colors crimson red navy blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup    --startLabel TSS  --endLabel TES  --plotHeight 9  --plotWidth 12 --yMin 0 --yMax 1  --yAxisLabel "meCpG level" -out WGBS_EMseq_gene_body_me_CpH.svg

#    plotProfile -m WGBS_EMseq_enhancer.scale.5k.me.test.gz  --colors crimson red navy blue --samplesLabel WGBS_Rep1 WGBS_Rep2 EM_Rep1 EM_Rep2 --perGroup  --startLabel enhancer --endLabel enhancer -out WGBS_EMseq_enhancers_me_CpG.svg --plotHeight 9  --plotWidth 9 --yMin 0 --yMax 0.83  --yAxisLabel "meC level" &
# plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar' --zMin -3 --zMax 3 --kmeans 4
#  plotHeatmap -m  WGBS_EMseq_Str_EnHr.rp.gz   -out WGBS_EMseq_EnHr.pdf --colorMap RdBu --whatToShow 'heatmap and colorbar'  --kmeans 4
