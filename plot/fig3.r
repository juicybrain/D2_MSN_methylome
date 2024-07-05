setwd("C:/Users/Yx Li/Documents/D2_methylome_project/WGBS_EM_comparison/")
set.seed(2)

library(tidyverse)


# mCpG level
mCpG_level <- c(0.805653,	0.807307,	0.80313,	0.811942)
sample=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2")
sample=factor(sample,levels=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2"))
method=c("AM-seq","AM-seq","EM-seq","EM-seq")
dat_meCpGLevel<- data.frame(sample, mCpG_level,method)
p1.1 <- ggplot(dat_meCpGLevel,aes(x=sample,y=mCpG_level,col=method)) + geom_col(width=0.4,size=1) + coord_cartesian(ylim = c(0.75,0.85),expand=F)+scale_color_manual(values = c("EM-seq"="blue","AM-seq"="red"))+xlab("")+ylab("mCpG level (%)")+expand_limits(x =c(0.5,4.5))

p1.2=p1.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(1,0.5,2,0.5), "cm"),
                                                 legend.title= element_blank(),
                                                 axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
                                                 legend.position =c(0.2,0.8),panel.grid.major=element_line(colour=NA),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())




#+ theme_classic(base_size = 20)+ theme(plot.margin = unit(c(2,2,2,2), "cm"),axis.text.x =element_text(size = 16,angle = 30,vjust=1,hjust=1),legend.position = "top",panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())

# mCpH level
dat_mCpH <- read.table("mCpH.level.f.txt",sep="\t",header=T)
#dat_mCpH$CpN <- factor(dat_mCpH$CpN,levels=c("CpH","CpA","CpT","CpC"))
 p2.1 <- ggplot(dat_mCpH,aes(x=ID,y=mCpH_level*100, col=method, fill=CpN))+geom_col(position=position_dodge(0.8),size=1.2,width=0.6)+ylab(" meCH level (%)") +scale_y_continuous(expand = c(0, 0),limits = c(0,4))+scale_fill_manual(values = c("mCpC"="grey35" ,"mCpA"="grey15","mCpH"="grey95","mCpT"="grey65"))+scale_color_manual(values = c("EM-seq"="blue","AM-seq"="red"))+xlab(label = "") 
  
  p2.2 <- p2.1 +
    theme_classic(base_size = 16)+ 
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
          legend.title= element_blank(),
          legend.position =c(0.8,0.7),
          legend.box = "vertical",
          axis.text.x =element_text(size = 12,angle = 45,vjust=1,hjust = 1),
          panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = NA,size=1,colour = "black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          #panel.background = element_rect(fill = "transparent",colour = NA),
          #plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank()) + scale_x_continuous(breaks=dat_mCpH$ID,labels =dat_mCpH$Group,expand = c(0, 0), limits = c(0,17))
    
    
    
    #theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
     #                                            legend.title= element_blank(),
     #                                            axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
     #                                            legend.position =c(0.2,0.8),panel.grid.major=element_line(colour=NA),
      #                                           panel.background = element_rect(fill = NA,size=1,colour = "black"),
     #                                            plot.background = element_rect(fill = "transparent",colour = NA),
      #                                           #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
        #                                         panel.grid.minor = element_blank())
  
  
  
  
  
  
 # +theme_classic(base_size = 20)+ theme(plot.margin = unit(c(1,1,1,1), "cm"),axis.text.x =element_text(size = 16,angle = 45,vjust=1,hjust = 1),legend.position ="top",panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank()) + scale_x_continuous(breaks=dat_mCpH$ID,labels =dat_mCpH$Group,expand = c(0, 0), limits = c(0,17))
  
# mCG mCH correlation
  dat_mCG_mCH_10kbin <- read.table("CG_CH_melevel.10kbin.txt")
  dat_mCG_mCH_cor <- cor(dat_mCG_mCH_10kbin[1:8],method="pearson")
  row.names(dat_mCG_mCH_cor) <- c("AM_Rep1_CpG","AM_Rep2_CpG","EM_Rep1_CpG","EM_Rep2_CpG","AM_Rep1_CpH","AM_Rep2_CpH","EM_Rep1_CpH","EM_Rep2_CpH")
  colnames(dat_mCG_mCH_cor) <- c("AM_Rep1_CpG","AM_Rep2_CpG","EM_Rep1_CpG","EM_Rep2_CpG","AM_Rep1_CpH","AM_Rep2_CpH","EM_Rep1_CpH","EM_Rep2_CpH")

  library(reshape2)
dat_CN_cor_ggplot <- melt(dat_mCG_mCH_cor)

p3.1 <-ggplot(dat_CN_cor_ggplot, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  geom_text(size=4,aes(label =round(value,2))) +
  scale_fill_gradientn(colours = c("green","blue","yellow","orange","red"),limits=c(0,1))+
  labs(x = "", y = "", title = "mCpG mCpH level correlation")+xlab(label = "")+ylab(label = "")+
  labs(fill= "r value")+
  theme_classic(base_size = 16)+ theme(plot.margin = unit(c(1.5,0.1,1.5,0.1), "cm"),panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank(),axis.text.x = element_text(vjust=1,hjust=1,angle = 45,size=12),axis.text.y = element_text(angle =0,size=12))



#top 10000 CH site
library(ggseqlogo)
 EM_Rep1_CpN <- read.table("EM_D2_Rep1.all.top10000.fa.txt")
p_motif_1 <- ggseqlogo(as.character(unlist(EM_Rep1_CpN$V1)),method="b") +theme(panel.background = element_rect(fill = NA,size=1,colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(1,0,0,0), "cm"))+xlab(label = "EM_Rep1")

 EM_Rep2_CpN <- read.table("EM_D2_Rep2.all.top10000.fa.txt")
 p_motif_2 <- ggseqlogo(as.character(unlist(EM_Rep2_CpN$V1)),method="b")+theme(panel.background = element_rect(fill = NA,size=1,colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(1,0,0,0), "cm"))+xlab(label = "EM_Rep2")
 AM_Rep1_CpN <- read.table("AM_D2_Rep1.all.top10000.fa.txt")
p_motif_3 <- ggseqlogo(as.character(unlist(AM_Rep1_CpN$V1)),method="b")+theme(panel.background = element_rect(fill = NA,size=1,colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(1,0,0,0), "cm")) +xlab(label = "AM_Rep1")

 AM_Rep2_CpN <- read.table("AM_D2_Rep2.all.top10000.fa.txt")
 p_motif_4 <- ggseqlogo(as.character(unlist(AM_Rep2_CpN$V1)),method="b")+theme(panel.background = element_rect(fill = NA,size=1,colour = "black"),plot.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(1,0,0,0), "cm")) +xlab(label = "AM_Rep2")

 library(grid)
 p.CpT <- gridExtra::grid.arrange(p_motif_1, p_motif_2,p_motif_3, p_motif_4,vp=viewport(width = 0.8, height = 0.8),layout_matrix = matrix(c(3,4,1,2), nrow =2))
 #+ theme(plot.margin = unit(c(1,0.5,0.5,1), "cm"))
 
library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

#figure <- ggarrange(plotlist =list(p2,p3,p4),labels=c("A","B","C"),ncol=2,nrow=1)

# figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 2, 3, 4,5,6), nrow = 6))
figure1 <- grid.arrange( p1.2, p2.2,p3.1,p.CpT,layout_matrix = matrix(c(1,3,2,NA,4,NA), nrow =2))

ggsave("figre3.test1.svg",plot=figure1, device="svg",width=44, height=30,units="cm")
