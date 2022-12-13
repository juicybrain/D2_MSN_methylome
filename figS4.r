library("ChIPseeker")
library("tidyverse")
library("clusterProfiler") 
library("org.Mm.eg.db") 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

dat_CGI_inside_UMR <- readPeakFile("C:/Users/Yx Li/Documents/D2_methylome_project/wgbs_EM_comparison/mm10.CGI.mappable_in_UMR.bed")
dat_CGI_outside_UMR <- readPeakFile("mm10.CGI.mappable_outside_of_UMR.bed")
dat_UMR_without_CGI <- readPeakFile("D2.UMR.without.CGI.bed")
dat_UMR_without_CGI_intergenic <- readPeakFile("D2.UMR.without.CGI.intergenic.bed")

peaks <- list(CGIinUMR=dat_CGI_inside_UMR, CGIoutUMR=dat_CGI_outside_UMR, UMRwithoutCGI=dat_UMR_without_CGI, UMRwithoutCGI_intergenic =dat_UMR_without_CGI_intergenic )

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 100), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Mm.eg.db")

dat_gene_list = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

GO_1 <- enrichGO(gene = dat_gene_list$CGIinUMR, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
GO_2 <- enrichGO(gene = dat_gene_list$CGIoutUMR, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.001, readable = TRUE)
GO_2_simp <-  clusterProfiler::simplify(GO_2, cutoff = 0.4, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

GO_3 <- enrichGO(gene = dat_gene_list$UMRwithoutCGI, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.001, readable = TRUE)
GO_3_simp <- clusterProfiler::simplify(GO_3, cutoff = 0.4, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

GO_4 <- enrichGO(gene = dat_gene_list$UMRwithoutCGI_intergenic, keyType = "ENTREZID", OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.001, readable = TRUE)
GO_4_simp <- clusterProfiler::simplify(GO_4, cutoff = 0.4, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

#GO_1_simp <- clusterProfiler::simplify(GO_1, cutoff = 0.7, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

GO_2_simp <- clusterProfiler::simplify(GO_2, cutoff = 0.7, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

compGO_BP <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "BP",pvalueCutoff=0.001,qvalueCutoff=0.001)
compGO_CC <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "CC",pvalueCutoff=0.001,qvalueCutoff=0.001)
compGO_MF <- compareCluster(dat_gene_list,fun="enrichGO",OrgDb = org.Mm.eg.db, ont = "MF",pvalueCutoff=0.001,qvalueCutoff=0.001)

compGO_simp_BP <- clusterProfiler::simplify(compGO_BP, cutoff = 0.5, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)
compGO_simp_CC <- clusterProfiler::simplify(compGO_CC, cutoff = 0.5, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)
compGO_simp_MF <- clusterProfiler::simplify(compGO_MF, cutoff = 0.5, by = "p.adjust",select_fun = min, measure = "Wang", semData = NULL)

#p5.1 <- dotplot(compGO_simp,showCategory = 50, title = "GO Analysis")
p6.1 <- dotplot(compGO_simp_BP,showCategory = 20, title = "GO Analysis")
p7.1 <- dotplot(compGO_simp_CC,showCategory = 10, title = "GO Analysis")
p8.1 <- dotplot(compGO_simp_MF,showCategory = 10, title = "GO Analysis")
ggsave("UMR.GO_BP.svg",plot=p6.1, device="svg",width=12, height=10)
ggsave("UMR.GO_CC.svg",plot=p7.1, device="svg",width=12, height=8)
ggsave("UMR.GO_MF.svg",plot=p8.1, device="svg",width=12, height=8)


p1.1 <- dotplot(GO_2_simp, showCategory = 50, title = "GO Analysis")
p2.1 <- dotplot(GO_3_simp,showCategory = 50, title = "GO Analysis")
p3.1 <- dotplot(GO_4_simp,showCategory = 50, title = "GO Analysis")

compKEGG <- compareCluster(geneCluster = dat_gene_list, fun = "enrichKEGG", organism = "mouse", pvalueCutoff = 0.01, pAdjustMethod = "BH")
p4.1 <- dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
ggsave("UMR.KEGG.svg",plot=p4.1, device="svg",width=12, height=8)

ComB <- function(x){paste(bitr(unlist(strsplit(x,"/")), fromType = "ENTREZID",toType = c("ENSEMBL","SYMBOL"), OrgDb = org.Mm.eg.db)$SYMBOL,collapse = "/")}
GO_cluster_summary <- as.data.frame(compGO_simp_BP)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "GO_UMR_BP.txt",sep="\t",quote=F,row.names = F)

GO_cluster_summary <- as.data.frame(compGO_simp_CC)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "GO_UMR_CC.txt",sep="\t",quote=F,row.names = F)

GO_cluster_summary <- as.data.frame(compGO_simp_MF)
GO_cluster_summary$Symbol <- sapply(GO_cluster_summary$geneID,FUN=ComB)
write.table(GO_cluster_summary, "GO_UMR_MF.txt",sep="\t",quote=F,row.names = F)

KEGG_cluster_summary <- as.data.frame(compKEGG)
KEGG_cluster_summary$Symbol <- sapply(KEGG_cluster_summary$geneID,FUN=ComB)
write.table(KEGG_cluster_summary, "KEGG_UMR.txt",sep="\t",quote=F,row.names = F)


KEGG_4 <- enrichKEGG(gene = dat_gene_list$UMRwithoutCGI_intergenic, organism = "mouse" , pAdjustMethod = "BH", qvalueCutoff = 0.01)

p5.1 <- dotplot(KEGG_4, showCategory = 50, title = "GO Analysis")
ggsave("UMR_intergenic_region.KEGG.svg",plot=p5.1, device="svg",width=6, height=5)

KEGG_4_summary <- as.data.frame(KEGG_4)
KEGG_4_summary$Symbol <- sapply(KEGG_4_summary$geneID,FUN=ComB)
write.table(KEGG_4_summary, "KEGG_UMR_intergenic_region.txt",sep="\t",quote=F,row.names = F)
