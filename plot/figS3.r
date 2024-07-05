setwd("C:/Users/Yx Li/Documents/D2_methylome_project/wgbs_EM_comparison/")
set.seed(2)

library(tidyverse)
 library(grid)
 library(ggseqlogo)
library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

# CpT   motif

library(ggseqlogo)
 EM_Rep1_CpT <- read.table("EM_D2_Rep1.CpT.top10k.fa")
p_motif_1 <- ggseqlogo(as.character(unlist(EM_Rep1_CpT$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 EM_Rep2_CpT <- read.table("EM_D2_Rep2.CpT.top10k.fa")
 p_motif_2 <- ggseqlogo(as.character(unlist(EM_Rep2_CpT$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))
 wgbs_Rep1_CpT <- read.table("wgbs_D2_Rep1.CpT.top10k.fa")
p_motif_3 <- ggseqlogo(as.character(unlist(EM_Rep1_CpT$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 wgbs_Rep2_CpT <- read.table("wgbs_D2_Rep2.CpT.top10k.fa")
 p_motif_4 <- ggseqlogo(as.character(unlist(EM_Rep2_CpT$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 p.CpT <- gridExtra::grid.arrange(p_motif_1, p_motif_2,p_motif_3, p_motif_4,vp=viewport(width = 0.9, height = 0.9),layout_matrix = matrix(c(3,4,1,2), nrow =2))


library(ggseqlogo)
 EM_Rep1_CpA <- read.table("EM_D2_Rep1.CpA.cov.gz.top10000.fa.txt")
p_motif_1 <- ggseqlogo(as.character(unlist(EM_Rep1_CpA$V1)),method="b") +theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 EM_Rep2_CpA <- read.table("EM_D2_Rep2.CpA.cov.gz.top10000.fa.txt")
 p_motif_2 <- ggseqlogo(as.character(unlist(EM_Rep2_CpA$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))
 wgbs_Rep1_CpA <- read.table("wgbs_D2_Rep1.CpA.cov.gz.top10000.fa.txt")
p_motif_3 <- ggseqlogo(as.character(unlist(wgbs_Rep1_CpA$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 wgbs_Rep2_CpA <- read.table("wgbs_D2_Rep2.CpA.cov.gz.top10000.fa.txt")
 p_motif_4 <- ggseqlogo(as.character(unlist(wgbs_Rep2_CpA$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 library(grid)
 p.CpA <- gridExtra::grid.arrange(p_motif_1, p_motif_2,p_motif_3, p_motif_4,vp=viewport(width = 0.9, height = 0.9),layout_matrix = matrix(c(3,4,1,2), nrow =2))
 
 
 library(ggseqlogo)
 EM_Rep1_CpC <- read.table("EM_D2_Rep1.CpC.GT7.top1000.fa.txt")
p_motif_1 <- ggseqlogo(as.character(unlist(EM_Rep1_CpC$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 EM_Rep2_CpC <- read.table("EM_D2_Rep2.CpC.GT7.top1000.fa.txt")
 p_motif_2 <- ggseqlogo(as.character(unlist(EM_Rep2_CpC$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))
 wgbs_Rep1_CpC <- read.table("wgbs_D2_Rep1.CpC.GT7.top1000.fa.txt")
p_motif_3 <- ggseqlogo(as.character(unlist(wgbs_Rep1_CpC$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 wgbs_Rep2_CpC <- read.table("wgbs_D2_Rep2.CpC.GT7.top1000.fa.txt")
 p_motif_4 <- ggseqlogo(as.character(unlist(wgbs_Rep2_CpC$V1)),method="b")+theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm"))

 library(grid)
 p.CpC <- gridExtra::grid.arrange(p_motif_1, p_motif_2,p_motif_3, p_motif_4,vp=viewport(width = 0.9, height = 0.9),layout_matrix = matrix(c(3,4,1,2), nrow =2))
 
 p1.1 <- grid.arrange( p.CpA, p.CpT,p.CpC,layout_matrix = matrix(c(1,2,3), nrow =1),vp=viewport(width = 0.9, height = 0.9))

ggsave("Sfigre3.test1.svg",plot=figure1, device="svg",width=24, height=6)
 
 
# meCpN composition


cx_dat <- read.table("CpN.txt",sep="\t",header=T)
cx_dat$meCpN <- factor(cx_dat$meCpN,levels=c("meCG","meCA","meCC","meCT"))



library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 4))

 p2.1 <- ggplot(cx_dat,aes(x=sample,fill=meCpN,y=percentage)) + geom_bar(stat='identity',width = 0.5)+coord_flip()+scale_fill_manual(values =c( "meCG"="#1B9E77", "meCA"="#D95F02", "meCC"="#7570B3" ,"meCT"="#E7298A"))+xlab("")+ylab("meCpN/meC (%)") +geom_text(size=8,aes(sample,label=round(percentage,1)))
 
p2.2 <- p2.1 +theme_classic(base_size = 20)+ theme(plot.margin = unit(c(2,2,2,2), "cm"),axis.text.x =element_text(size = 16,angle = 45,vjust=0.5),legend.position ="top",panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())

# meCpN sites compostion


cx_dat_binom <- read.table("CpN_binom.txt",sep="\t",header=T)
cx_dat_binom$meCpN_sites <- factor(cx_dat_binom$meCpN_sites,levels=c("meCG_sites","meCA_sites","meCC_sites","meCT_sites"))



library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 4))

 p3.1 <- ggplot(cx_dat_binom,aes(x=sample,fill=meCpN_sites,y=percentage)) + geom_bar(stat='identity',width = 0.5)+coord_flip()+scale_fill_manual(values =c( "meCG_sites"="#1B9E77", "meCA_sites"="#D95F02", "meCC_sites"="#7570B3" ,"meCT_sites"="#E7298A"))+xlab("")+ylab("meCpN/meC (%)") +geom_text(size=8,aes(sample,label=round(percentage,2)))
 
p3.2 <- p3.1 +theme_classic(base_size = 20)+ theme(plot.margin = unit(c(2,2,2,2), "cm"),axis.text.x =element_text(size = 16,angle = 45,vjust=0.5),legend.position ="top",panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())





figure1 <- grid.arrange( p1.1, p2.2,p3.2,layout_matrix = matrix(c(2,1,3,1), nrow =2))

ggsave("Sfigre3.test1.svg",plot=figure1, device="svg",width=24, height=12)
