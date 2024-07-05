setwd("C:/Users/Yx Li/Documents/D2_methylome_project/WGBS_EM_comparison/")
set.seed(2)

library(tidyverse)

#CpG coverage fraction

dat_CpG_covRate <- read.table("AM_Rep1.CpG.cov.rate")
dat_CpG_covRate <- dat_CpG_covRate[,-1]
colnames <- c("CpG_coverage","Percentage")
dat_CpG_covRate$sample <- "AM_Rep1"

for (i in c("AM_Rep2","EM_Rep1","EM_Rep2")){

  path <- paste0(i,".CpG.cov.rate")
  dat_tmp <- read.table(path,header=F)
  dat_tmp <- dat_tmp[,-1]
  dat_tmp$sample <- i
  dat_CpG_covRate <- rbind(dat_CpG_covRate, dat_tmp)
}

p1.1 <- ggplot(dat_CpG_covRate,aes(x=V2,y=100*V3,col=sample))+ geom_line(size=1)+xlab("CpG coverage")+ylab("CpG over x fold coverage(%)")+scale_color_manual(values=c("EM_Rep1"="blue4","EM_Rep2"="blue","AM_Rep1"="red","AM_Rep2"="violetred"))+scale_y_continuous(expand = c(0, 0),limits=c(0,100)) +scale_x_continuous(expand = c(0, 0),limits=c(1,15),labels = as.character(dat_CpG_covRate$V2), breaks = dat_CpG_covRate$V2)

p1.2 <- p1.1 +theme_classic(base_size = 16)+ 
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        legend.title= element_blank(),
        axis.text.x =element_text(size = 12),
        legend.position =c(0.8,0.8),panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = NA,size=1,colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        #panel.background = element_rect(fill = "transparent",colour = NA),
        #plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

  #theme_classic(base_size = 24)+ theme(plot.margin = unit(c(2,2,2,2), "cm"),axis.text.x =element_text(size = 18),legend.position = c(0.8, 0.6),panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())
  
  
  #theme_classic(base_size = 26)+ theme(plot.margin = unit(c(1,1,1,1), "cm"),axis.text.x =element_text(size = 24),legend.position ="top",legend.title =element_text(color = "black", size = 14),legend.text = element_text(color = "black", size =12),panel.grid.major=element_line(colour=NA),panel.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA),panel.grid.minor = element_blank())

# coverage per dinucleotides

#dat <- read.table("dinucleotide_test.txt",header=T)

#p1 <- ggplot(dat,aes(x=,y=va))+ geom_line(aes(color=Group),size=1)+xlab("GC")+ylab("normalized coverage")

#p1=p1 + theme_grey(base_size = 12)+theme(legend.position = "top")

#ggsave("GC_content.pdf",p1,width = 15,height = 7)

#dat$x <- row.names(dat)
#pdf("dinucleotide.pdf",width=10,height = 7)
#par(pin=c(9,6))
#par(lwd=2,cex=0.5)
#par(cex.axis=1.5, cex.lab=3,font.axis=1,font.lab=1)
#plot(x=dat$x,y=dat$AM_Rep1,type ="l",xaxt = "n",col="red3",ylim=c(-0.01,0.01),lwd=2,xlab = "dinucleotide",ylab = "variation")
#lines(x=dat$x,y=dat$AM_Rep2, col = "red",lwd=2,pch=18)
#lines(x=dat$x,y=dat$EM_Rep1, col = "navy",lwd=2,pch=18)
##lines(x=dat$x,y=dat$EM_Rep2, col = "blue",lwd=2,pch=18)
#legend(1,0.3,legend=c("EMseq_Rep1","EMseq_Rep2","EM_Rep1","EM_Rep2"),col=c("navy", "blue","red3","red"),lty=1:2,cex=1.5)
#axis(1,at=1:16,labels=dat$Dinucleotide)
#dev.off()

dat_diyad_cov <- read.table("dinucleotide_test.txt",sep="\t", header=T)
dat_diyad_cov$Dinucleotide <- as.factor(dat_diyad_cov$Dinucleotide) #,level=c("AA","AC","AG","AT","CA","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"))
p2.1 <- ggplot()+geom_hline(yintercept=0, size=1,linetype="dashed",alpha=12, color = "black")+ geom_line(data=dat_diyad_cov,aes(x=group,y=coverage_rate_deviation,color=sample),size=1)+ylab("Coverage deviation(obs/bgd)")+ scale_color_manual(values=c("EM_Rep1"="blue4","EM_Rep2"="blue","AM_Rep1"="red","AM_Rep2"="violetred"))+scale_y_continuous(expand = c(0, 0),limits=c(-0.2,0.2))+xlab(label="Dinucleotide")

x_lab <- read.table("diyad_x_label.txt")

p2.2=p2.1 +theme_classic(base_size = 16)+
    theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
          legend.title= element_blank(),
          axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
          legend.position = c(0.8,0.75),panel.grid.major=element_line(colour=NA),
          #panel.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = NA,size=1,colour = "black"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          #plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())+scale_x_continuous(breaks=x_lab$V2,
          labels = x_lab$V1,expand = c(0, 0), limits = c(0,17))
  
  
  
  
  #theme(plot.margin = unit(c(1,1,1,1), "cm"),
  #axis.text.x =element_text(hjust=-0.4,vjust=-0.1,angle=45,size = 16),
  #legend.text=element_text(size=12),
  #legend.position = c(0.8, 0.8),panel.grid.major=element_line(colour=NA),
  #panel.background = element_rect(fill = NA,size=2,colour = "black"),
  #plot.background = element_rect(fill = "transparent",colour = NA),
  #panel.grid.minor = element_blank())+
  #scale_x_continuous(breaks=x_lab$V2,labels = x_lab$V1,expand = c(0, 0), limits = c(0,17))


# coverge per GC content

dat_GC <- read.table("EM_Rep1.gc_bias.txt",header=F,sep="\t",skip=7,comment.char = "#")

dat_GC$sample <- "EM_Rep1"
for (i in c("EM_Rep2","AM_Rep1","AM_Rep2")){

  path <- paste0(i,".gc_bias.txt")
  dat_tmp <- read.table(path,header=F,sep="\t",skip=7,comment.char = "#")
  dat_tmp$sample <- i
 # dat_tmp$V2 <- 1e8*dat_tmp$V2/sum(dat_tmp$V2)
  dat_GC <- rbind(dat_GC, dat_tmp)
  }
#bin sample into 5 bp binsize
df_trans <- transform(dat_GC,group=cut(V3,breaks=c(seq(-1,100,5)),labels=c(seq(5,100,5))))
#df_trans[is.na(df_new)] <- 0
df_trans$V12 <- paste(df_trans$sample,df_trans$group,sep="-")
df_res <- do.call(data.frame,aggregate(V7~V12,df_trans,FUN=function(x) c(Count=length(x),Sum=sum(x))))
df_new <- separate(data = df_res, col = V12, into = c("sample", "bp"), sep = "-")
df_new[is.na(df_new)] <- 0
df_new$size <- as.numeric(df_new$bp)-2.5

p3.1 <- ggplot(df_new,aes(x=size,y=V7.Sum))+ geom_line(aes(color=sample),size=1)+xlab("CG content(%)")+ylab("Normalized coverage")+ scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 10))+scale_color_manual(values=c("EM_Rep1"="blue4","EM_Rep2"="blue","AM_Rep1"="red","AM_Rep2"="violetred"))

p3.2=p3.1+theme_classic(base_size = 16)+ theme(plot.margin = unit(c(1,1,0.5,0.5), "cm"),
                                               legend.title= element_blank(),
                                               axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
                                               legend.position = c(0.8,0.7),panel.grid.major=element_line(colour=NA),
                                               #panel.background = element_rect(fill = "transparent",colour = NA),
                                               #plot.background = element_rect(fill = "transparent",colour = NA),
                                               panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                               plot.background = element_rect(fill = "transparent",colour = NA),
                                               panel.grid.minor = element_blank())

# low CC coverage in AM libraries
  dat_CC_cov <- read.table("CpN coverage.txt",header = T)
  
  p4.1 <- ggplot(dat_CC_cov,aes(x=ID,y=mean_coverage, col=method, fill=CpN))+geom_col(position=position_dodge(0.8),size=1.2,width=0.6)+ylab(" Mean coverage") +scale_y_continuous(expand = c(0, 0))+scale_fill_manual(values = c("CpA"="grey15","CpC"="grey95","CpG"="grey35" ,"CpT"="grey65"))+scale_color_manual(values = c("EM-seq"="blue","AM"="red"))+xlab(label = "") + coord_cartesian(ylim=c(5,10))
  
  p4.2 <- p4.1 +theme_classic(base_size = 16)+ theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),legend.title= element_blank(),
                                                     axis.text.x =element_text(size = 12,angle = 45,vjust=0.5),
                                                     legend.position =c(0.2,0.8),
                                                     legend.box = "horizontal",
                                                     panel.grid.major=element_line(colour=NA),
                                                     panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                     plot.background = element_rect(fill = "transparent",colour = NA),
                                                     #panel.background = element_rect(fill = "transparent",colour = NA),
                                                     #plot.background = element_rect(fill = "transparent",colour = NA),
                                                     panel.grid.minor = element_blank()) +
                                                     scale_x_continuous(breaks=dat_CC_cov$ID,labels =dat_CC_cov$Group,expand = c(0, 0), limits = c(0,17))
  
  
  
  
# combined together
library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

#figure <- ggarrange(plotlist =list(p2,p3,p4),labels=c("A","B","C"),ncol=2,nrow=1)

# figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 2, 3, 4,5,6), nrow = 6))
figure1 <- grid.arrange( p1.2, p2.2,p3.2,p4.2,layout_matrix = matrix(c(1,NA,2,4,3,4), nrow =2))

ggsave("figre2.test2.svg",plot=figure1, device="svg",width=44, height=30,unit="cm")

