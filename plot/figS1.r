setwd("C:/Users/Yx Li/Documents/D2_methylome_project/wgbs_EM_comparison/")
set.seed(1) 
library(tidyverse)
library(ggseqlogo)

# wgbs_D2_2_S20_R1_val_1_bismark_bt2_pe.deduplicated.cytosine_context_summary.txt unmethylated motif

 unconvtC <- read.table("wgbs_D2_2_S20_R1_val_1_bismark_bt2_pe.deduplicated.cytosine_context_summary.txt")  
 unconvtSeq <- c(rep(unconvtC$V1,unconvtC$V2))
 UNCONVTSEQ_sorted_w2 <- unconvtC[sort(unconvtC$V2,index.return=TRUE,decreasing = T)$ix,]
 dat_TMP <- c()
 dat_motif <- c()
 for (i in 1:64 ){
    dat_TMP <- c(rep(as.character(unlist(UNCONVTSEQ_sorted_w2[i,][1])), as.numeric(unlist(UNCONVTSEQ_sorted_w2[i,][2]))))
    dat_motif <- c(dat_motif, dat_TMP)
    
 }
 dat_motif_w2 <- dat_motif
 p_motif_1 <- ggseqlogo(dat_motif_w2,method="p")+theme(plot.margin = unit(c(1,1,1,1), "cm"))
 
 unconvtC <- read.table("wgbs_D2_1_S19_R1_val_1_bismark_bt2_pe.deduplicated.cytosine_context.s.txt")  
 unconvtSeq <- c(rep(unconvtC$V1,unconvtC$V2))
 UNCONVTSEQ_sorted_w1 <- unconvtC[sort(unconvtC$V2,index.return=TRUE,decreasing = T)$ix,]
 dat_TMP <- c()
 dat_motif <- c()
 for (i in 1:64 ){
    dat_TMP <- c(rep(as.character(unlist(UNCONVTSEQ_sorted_w1[i,][1])), as.numeric(unlist(UNCONVTSEQ_sorted_w1[i,][2]))))
    dat_motif <- c(dat_motif, dat_TMP)
    
 }
 dat_motif_w1 <- dat_motif
 
  p_motif_2 <- ggseqlogo(dat_motif_w1,method="p")+theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  
   unconvtC <- read.table("EM_D2_1_S17_R1_val_1_bismark_bt2_pe.deduplicated.cytosine_context.txt")  
 unconvtSeq <- c(rep(unconvtC$V1,unconvtC$V2))
 UNCONVTSEQ_sorted_E1 <- unconvtC[sort(unconvtC$V2,index.return=TRUE,decreasing = T)$ix,]
 dat_TMP <- c()
 dat_motif <- c()
 for (i in 1:64 ){
    dat_TMP <- c(rep(as.character(unlist(UNCONVTSEQ_sorted_E1[i,][1])), as.numeric(unlist(UNCONVTSEQ_sorted_E1[i,][2]))))
    dat_motif <- c(dat_motif, dat_TMP)
    
 }
 dat_motif_E1 <- dat_motif
 p_motif_3 <- ggseqlogo(dat_motif_E1,method="p")+theme(plot.margin = unit(c(1,1,1,1), "cm"))
 
  
 unconvtC <- read.table("EM_D2_2_S18_R1_val_1_bismark_bt2_pe.deduplicated.cytosine_context.txt")  
 unconvtSeq <- c(rep(unconvtC$V1,unconvtC$V2))
 UNCONVTSEQ_sorted_E2 <- unconvtC[sort(unconvtC$V2,index.return=TRUE,decreasing = T)$ix,]
 UNCONVTSEQ_sorted_E2$V2 <- UNCONVTSEQ_sorted_E2$V2*100
 dat_TMP <- c()
 dat_motif <- c()
 for (i in 1:64){
    dat_TMP <- c(rep(as.character(unlist(UNCONVTSEQ_sorted_E2[i,][1])), as.numeric(unlist(UNCONVTSEQ_sorted_E2[i,][2]))))
    dat_motif <- c(dat_motif, dat_TMP)
    
 }
 dat_motif_E2 <- dat_motif
 p_motif_4 <- ggseqlogo(dat_motif_E2,method="p")+theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
 p_unconver_all <-  gridExtra::grid.arrange(p_motif_1, p_motif_2,p_motif_3,p_motif_4,vp=grid::viewport(width=0.9,height=0.7))


 # duplication rate correlation
 
 sample=c("WGBS_Rep1","WGBS_Rep2","EM_Rep1","EM_Rep2")
sample=factor(sample,levels=c("WGBS_Rep1","WGBS_Rep2","EM_Rep1","EM_Rep2"))
rate=c(0.1695,0.1747,0.1351,0.1505)
method=c("WGBS","WGBS","EM-seq","EM-seq")
survivalRate <- data.frame(sample, rate,method)
p1.1 <- ggplot(survivalRate,aes(x=sample,y=rate,fill=method)) + geom_col(width=0.4)+scale_fill_manual(values = c("WGBS"="grey40","EM-seq"="grey65"))+xlab("samples")+ylab("Duplication rate")+scale_y_continuous(expand = c(0, 0),limits=c(0,0.3))
p_dup=p1.1 + theme_classic(base_size = 24)+ theme(plot.margin = unit(c(2,2,2,2), "cm"),axis.text.x =element_text(size = 18,angle = 30,vjust=0.5),legend.position = "top",panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())
 

 
 
 # coverage blocks correlation
 library("pheatmap")
 library("RColorBrewer")
 
 
dat_cvg_cor <- read.table("coverage_cor.txt")

dat_cor <- cor(dat_cvg_cor[2:5],method="spearman")
row.names(dat_cor) <- c("EM_Rep1","EM_Rep2","WGBS_Rep1","WGBS_Rep2")
colnames(dat_cor) <- c("EM_Rep1","EM_Rep2","WGBS_Rep1","WGBS_Rep2")



bk <- c(seq(0,1,by=0.01))
p3 <- pheatmap(dat_cor,display_numbers=T,fontsize_number = 8,legend_breaks = seq(0,1,0.2),breaks=bk,treeheight_row = 0, treeheight_col = 0,annotation_names_col=F)

#use ggplot2
library(reshape2)
dat_cor_ggplot <- melt(dat_cor)

p_cor <-ggplot(dat_cor_ggplot, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  geom_text(size=8,aes(label =round(value,2))) +
  scale_fill_gradientn(colours = c("white","yellow","blue"),limits=c(0,1))+
  labs(x = "", y = "", title = "coverage correlation(spearman)")+
  labs(fill= "rho value")+
  theme_classic(base_size = 24)+ theme(plot.margin = unit(c(1,0,1,1), "cm"),panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(),axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45,size=18))
                                  

 
library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

#figure <- ggarrange(plotlist =list(p2,p3,p4),labels=c("A","B","C"),ncol=2,nrow=1)

# figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 2, 3, 4,5,6), nrow = 6))
figure1 <- grid.arrange( p_unconver_all, p_dup,p_cor,layout_matrix = matrix(c(1, 2,1,3), nrow = 2))

ggsave("Sfigure1.svg",plot=figure1, device="svg",width=24, height=24)
    
 
 
