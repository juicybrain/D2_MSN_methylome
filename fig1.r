setwd("C:/Users/Yx Li/Documents/D2_methylome_project/WGBS_EM_comparison/")
set.seed(2)

library(tidyverse)

# reads survival rate
sample=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2")
sample=factor(sample,levels=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2"))
rate=c(0.9971,0.9973,0.9992,0.9993)
method=c("AM-seq","AM-seq","EM-seq","EM-seq")
survivalRate <- data.frame(sample, rate,method)
p1.1 <- ggplot(survivalRate,aes(x=sample,y=rate*100,col=method)) + geom_col(width=0.4,size=1) + 
  coord_cartesian(ylim = c(99,100),expand=F)+expand_limits(x =c(0.5,4.5))+
  #scale_y_continuous(limits = c(99, 100),expand=c(0,0))+
  scale_color_manual(values = c("AM-seq"="red","EM-seq"="blue"))+xlab("")+ylab("Reads survival rate(%)")
p1.2=p1.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
                                                 axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
                                                 legend.position ="none",panel.grid.major=element_line(colour=NA),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                 #plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

# conversion rate
#    \u03bb   lambda ??

#library("ggpattern")
dat_conversion <- read.table("WGBS_EMseq_conversion.txt",header=T)
dat_conversion$sample <- factor(dat_conversion$sample,levels=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2"))
p2.1 <- ggplot(dat_conversion,aes(x=sample,y=coversion_rate*100, col=method, fill=C_context))+geom_col(position=position_dodge(0.8),size=1,width=0.6)+xlab("")+ylab("Unconversted C of unmethylated \u03bb-DNA (%)") +scale_y_continuous(expand = c(0, 0),limits=c(0,0.4))+scale_fill_manual(values = c("CpA"="grey15","CpC"="grey95","CpG"="grey35" ,"CpT"="grey65"))+scale_color_manual(values = c("EM-seq"="blue","AM-seq"="red"))
 
p2.2=p2.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(1,1,0.5,0.5), "cm"),
                                                 axis.text.x =element_text(size = 12),
                                                 legend.position = c(0.8, 0.6),panel.grid.major=element_line(colour=NA),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                # panel.background = element_rect(fill = "transparent",colour = NA),
                                                # plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())
 
# meCpG conversion rate
sample=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2")
sample=factor(sample,levels=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2"))
meCpG_cvn_Rate <- c(0.0322795,0.0358482,0.00960995,0.0191819)
method <- c("AM-seq","AM-seq","EM-seq","EM-seq")
dat_meCpG_cvnRate <- data.frame(sample,meCpG_cvn_Rate,method)

p3.1 <- ggplot(survivalRate,aes(x=sample,y=meCpG_cvn_Rate*100, col=method)) + geom_col(width=0.4,size=1)+scale_color_manual(values = c("AM-seq"="red","EM-seq"="blue"))+xlab("")+ylab("Coversion rate of meCpG in PUC19 DNA(%) ")+scale_y_continuous(expand = c(0, 0),limits=c(0,4))
p3.2=p3.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                                                 axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
                                                 legend.position =c(0.8, 0.8),panel.grid.major=element_line(colour=NA),
                                                 panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                 plot.background = element_rect(fill = "transparent",colour = NA),
                                                # panel.background = element_rect(fill = "transparent",colour = NA),
                                                # plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

# mapping ratio

sample=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2")
sample=factor(sample,levels=c("AM_Rep1","AM_Rep2","EM_Rep1","EM_Rep2"))
mapping_Rate <- c(71.4, 77.9, 69.8, 76.4)
method <- c("AM-seq","AM-seq","EM-seq","EM-seq")
dat_mapping <- data.frame(sample,mapping_Rate,method)

p4.1 <- ggplot(dat_mapping,aes(x=sample,y=mapping_Rate,col=method)) + geom_col(width=0.4,size=1) +scale_color_manual(values = c("AM-seq"="red","EM-seq"="blue"))+xlab("")+ylab("Mapping rate to mm10 genome(%)") + scale_y_continuous(expand = c(0, 0),limits=c(0,100))
p4.2=p4.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                                                 axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),
                                                 legend.position = c(0.8,0.8),panel.grid.major=element_line(colour=NA),
                                                  panel.background = element_rect(fill = NA,size=1,colour = "black"),
                                                  plot.background = element_rect(fill = "transparent",colour = NA),
                                                 #panel.background = element_rect(fill = "transparent",colour = NA),
                                                # plot.background = element_rect(fill = "transparent",colour = NA),
                                                 panel.grid.minor = element_blank())

# coverage across genome

#dat_cov_acrGenome <- read.table(paste0("EM_Rep1","_coverage_across_reference.txt"),header=F)
 #dat_2$sample <- "Emseq_1"
 #dat_cov_acrGenome$sample <- "EM_Rep1"
 #dat_3$V2 <- 1e8*dat_2$V2/sum(dat_2$V2)
#for (i in c("EM_Rep2","AM_Rep1","AM_Rep2")){

  #path <- paste0(i,"_coverage_across_reference.txt")
  #dat_tmp <- read.table(path,header=F)
  #dat_tmp$sample <- i
  #dat_tmp$Group <- strsplit(as.character(i),split="_")[[1]][1]
  #dat_tmp$V2 <- 1e8*dat_tmp$V2/sum(dat_tmp$V2)
  #dat_cov_acrGenome <- rbind(dat_cov_acrGenome, dat_tmp)
  #}
dat_cov_acrGenome  <- read.table("whole_genome_coverage.txt",header=T)
#dat_cov_acrGenome$method <- str_split_fixed(dat_cov_acrGenome$sample,"\\_",2)[,1]
chr_lab <- read.table("chr_label.txt")
p5.1 <- ggplot(dat_cov_acrGenome,aes(x=loc,y=cov,group=sample,color=sample))+geom_line(size=1)+ylim(0,30)+ scale_y_continuous(expand = c(0, 0), limits = c(0, 30))+ylab("Average coverage")+xlab("")+scale_color_manual(values=c("EM_Rep1"="black","EM_Rep2"="blue4","AM_Rep1"="red","AM_Rep2"="pink"))

p5.2 <- p5.1 +theme_classic(base_size = 16)+
  theme(plot.margin = unit(c(0,1,0.5,0.5), "cm"),
  axis.text.x =element_text(hjust=-0.4,vjust=-0.5,angle=45,size = 12),
  legend.position = c(0.8, 0.8),panel.grid.major=element_line(colour=NA),
  panel.background = element_rect(fill = NA,size=1,colour = "black"),
  plot.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=chr_lab$V2,labels = chr_lab$V1,expand = c(0, 0), limits = c(0,2725537669))


# insert size
dat_2 <- read.table(paste0("EM_Rep1","_insert_size_histogram.txt"),header=F)
#dat_2$sample <- "Emseq_1"
dat_2$Group <- "EM_Rep1"
dat_2$V2 <- 1e8*dat_2$V2/sum(dat_2$V2)

for (i in c("EM_Rep2","AM_Rep1","AM_Rep2")){

  path <- paste0(i,"_insert_size_histogram.txt")
  dat_tmp <- read.table(path,header=F)
  dat_tmp$Group <- i
  #dat_tmp$Group <- strsplit(as.character(i),split="_")[[1]][1]
  dat_tmp$V2 <- 1e8*dat_tmp$V2/sum(dat_tmp$V2)
  dat_2 <- rbind(dat_2, dat_tmp)
  }
# bin reads by 20bp
df_trans <- transform(dat_2,group=cut(V1,breaks=c(seq(0,540,10)),labels=c(seq(10,540,10))))
df_trans$V5 <- paste(df_trans$Group,df_trans$group,sep="*")
df_res <- do.call(data.frame,aggregate(V2~V5,df_trans,FUN=function(x) c(Count=length(x),Sum=sum(x))))
df_new <- separate(data = df_res, col = V5, into = c("Group", "bp"), sep = "\\*")
df_new$size <- as.numeric(df_new$bp)
df_new$sample <- df_new$Group
df_new$method <- str_split_fixed(df_new$Group,"\\_",2)[,1]

p6.1 <- ggplot(df_new,aes(x=size,group=sample,y=V2.Sum,color=sample))+ geom_line(size=1)+xlab("")+ylab("Insertion count")+ scale_color_manual(values = c("black", "blue4","red","pink"))+scale_y_continuous(expand = c(0, 0),limits=c(0,6000000))

p6.2=p6.1 + theme_classic(base_size = 16)+ theme(plot.margin = unit(c(0,1,0.5,0.5), "cm"),
    axis.text.x =element_text(size = 12,angle = 0,vjust=0.5),legend.position = c(0.8, 0.8),
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = NA,size=1,colour = "black"),
    plot.background = element_rect(fill = "transparent",colour = NA),
   # panel.background = element_rect(fill = "transparent",colour = NA),
   # plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank())




# combined together

library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

#figure <- ggarrange(plotlist =list(p2,p3,p4),labels=c("A","B","C"),ncol=2,nrow=1)

# figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 2, 3, 4,5,6), nrow = 6))
figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 3, 5,2, 4,5,2,6,5), nrow = 3))

ggsave("figre1.test4.svg",plot=figure1, device="svg",width=44, height=44,units="cm")
