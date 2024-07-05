#!/home/yli/bin/Rscript3.6

library(methylKit)

        file.list = list("wgbs_1_CpG.txt","wgbs_2_CpG.txt","wgbs_5_CpG.txt", "wgbs_6_CpG.txt","EMseq_1_CpG.txt","EMseq_2_CpG.txt","EMseq_3_CpG.txt","EMseq_4_CpG.txt","EM_11_S23.pe.s.bam_CpG.txt","EM_12_S24.pe.s.bam_CpG.txt","EM_1_S13.pe.s.bam_CpG.txt","EM_2_S14.pe.s.bam_CpG.txt","EM_3_S15.pe.s.bam_CpG.txt","EM_4_S16.pe.s.bam_CpG.txt","EM_5_S17.pe.s.bam_CpG.txt","EM_6_S18.pe.s.bam_CpG.txt","EM_7_S19.pe.s.bam_CpG.txt","EM_8_S20.pe.s.bam_CpG.txt","EM_9_S21.pe.s.bam_CpG.txt")
myobj = methRead(file.list,
            sample.id=list("wgbs_1","wgbs_2","wgbs_5","wgbs_6", "EMseq_1", "EMseq_2","EMseq_3","EMseq_4","EM_11_L","EM_12_L","EM_1_S","EM_2_S","EM_3_H","EM_4_H","EM_5_L","EM_6_L","EM_7_S","EM_8_S","EM_9_H"),
            assembly="mm10",
           treatment=c(rep(1,8),rep(0,11)),
            context="CpG",
            dbtype="tabix")

#getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)

#getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)


meth=unite(myobj, destrand=FALSE)


save(meth,file="wgbs_emseq_meth_HLS.RData")
pdf("wgbs_emseq_correlation.pdf",width=10, height=10)
getCorrelation(meth,plot=TRUE)
dev.off()
pdf("sample_cluster_HLS.pdf",width=10,height=10)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
dev.off()
pdf("sample_PCA_HLS.pdf", width=10,height=10)
pca= PCASamples(meth,scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,sd.threshold=0.5,filterByQuantile=TRUE,obj.return=TRUE)
dev.off()
print(pca)
library(tidyverse)
summary(pca)
        pca_dat_B4 <- pca$x[,1:2]
        pca_dat_B4 <- as.data.frame(pca_dat_B4)
        pca_dat_B4$name = c("wgbs_1","wgbs_2","wgbs_5","wgbs_6", "EMseq_1", "EMseq_2","EMseq_3","EMseq_4","EM_11_L","EM_12_L","EM_1_S","EM_2_S","EM_3_H","EM_4_H","EM_5_L","EM_6_L","EM_7_S","EM_8_S","EM_9_H")

        pca_dat_B4$group <- c('wgbs','wgbs','wgbs','wgbs','EMseq','EMseq','EMseq','EMseq',"EM_L","EM_L","EM_S","EM_S","EM_H","EM_H","EM_L","EM_L","EM_S","EM_S","EM_H")
        print(pca_dat_B4)
        print(pca$sdev)

        print(pca$sdev[1]/sum(pca$sdev))
        print(pca$sdev[2]/sum(pca$sdev))


        write.table(pca_dat_B4, "wgbs_emseq_HLS_PCA.txt", quote=F)

        pdf("gglot.wgbs_EM_HLS_test.pdf", width=10, height=10)
        ggplot(pca_dat_B4, aes(x=pca_dat_B4$PC1,y=pca_dat_B4$PC2,label=pca_dat_B4$name)) + geom_point(aes(color=group)) + geom_text(aes(label=name),hjust=-0.5,vjust=-0.5,size=2) + xlab(paste0("PC1"," 30%"))+ylab(paste0("PC2", " 12%"))+xlim(-40,70)+ylim(-30,50)+theme_bw(base_size=2)+ theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())
        dev.off()

