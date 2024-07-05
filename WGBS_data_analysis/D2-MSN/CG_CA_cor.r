set.seed(2)
library(tidyverse)
setwd("C:/Users/Yx Li/Documents/D2_methylome_project/wgbs_EM_comparison")
dat_mCA <- read.table("meCA_DMR.delta.txt",header=F)
colnames(dat_mCA) <- c("chr","pos","meCA")

dat_mCG <- read.table("D2-Excit.CG.f.txt",header=F)
colnames(dat_mCG) <- c("chr","pos","meCG")

dat <- merge(dat_mCG,dat_mCA,by=c("chr","pos"),all=F)

lmFit <- lm(meCG~meCA,data=dat)
p=ggplot(dat,aes(x=meCG,y=meCA))+geom_point()

dat["Q1"] <- rank(dat$meCG)
dat["Q2"] <- rank(dat$meCA)

p=ggplot(dat,aes(x=Q1,y=Q2))+geom_point()
p
pairs(dat)
quantile(dat$meCG) %>% class()

dat_E1 <- read.table("E1.CG.CA.txt")
colnames(dat_E1) <- c("E1","mCG","mCA")
p1=ggplot(dat_E1,aes(x=mCG,y=mCA))+geom_point(alpha=0.2) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4))
lm_FIt <- lm(mCA~mCG,data=dat_E1)

lm_FIt_2 <- lm(mCA~mCG+I(mCG^2),data=dat_E1)
lm_FIt_3 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3),data=dat_E1)
lm_FIt_4 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_E1)
#lm.fit5=lm(mCA???poly(mCG ,5),data=dat_E1)
lm.fit6=lm(mCA~log(mCG+0.1),data=dat_E1)

dat_E1_LT0.5 <- dat_E1%>%filter(mCG<0.5)
p1.1=ggplot(dat_E1_LT0.5,aes(x=mCG,y=mCA))+geom_point(alpha=0.2)
cor(dat_E1_LT0.5$mCG,dat_E1_LT0.5$mCA)
dat_E1_GT0.5 <- dat_E1%>%filter(mCG>0.5)
p1.2=ggplot(dat_E1_GT0.5,aes(x=mCG,y=mCA))+geom_point(alpha=0.2)
cor(dat_E1_GT0.5$mCG,dat_E1_GT0.5$mCA)
p1=ggplot(dat_E1,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG level in Excitatory neuron")+ylab("D2 CpG-DMR mCA level in Excitatory neuron")
lm_FIt_E1 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_E1)
#lm_FIT_10_E1=lm(mCA???poly(mCG ,10),data=dat_E1)

dat_E2 <- read.table("E2.CG.CA.txt")
colnames(dat_E2) <- c("E2","mCG","mCA")
colnames(dat_E2) <- c("E2","mCG","mCA")
p2=ggplot(dat_E2,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG level in Excitatory neuron")+ylab("D2 CpG-DMR mCA level in Excitatory neuron")
lm_FIt_E2 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_E2)
#lm_FIT_10_E2=lm(mCA???poly(mCG ,10),data=dat_E2)


dat_EM_1 <- read.table("EM_D2_Rep1.CG.CA.txt")
colnames(dat_EM_1) <- c("EM1","mCG","mCA")
p3=ggplot(dat_EM_1,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG in D2-MSN(EM-seq)")+ylab("D2 CpG-DMR mCA level in D2-MSN(EM-seq)")
lm_FIt_EM_1 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_EM_1)
#lm_FIT_10_EM1=lm(mCA???poly(mCG ,10),data=dat_EM_1)

dat_EM_2 <- read.table("EM_D2_Rep2.CG.CA.txt")
colnames(dat_EM_2) <- c("EM1","mCG","mCA")
p4=ggplot(dat_EM_2,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG in D2-MSN(EM-seq)")+ylab("D2 CpG-DMR mCA level in D2-MSN(EM-seq)")
lm_FIt_EM_2 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_EM_2)
#lm_FIT_10_EM2=lm(mCA???poly(mCG ,10),data=dat_EM_2)

dat_WG_1 <- read.table("wgbs_D2_Rep1.CG.CA.txt")
colnames(dat_WG_1) <- c("EM1","mCG","mCA")
p5=ggplot(dat_WG_1,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG in D2-MSN(AM-seq)")+ylab("D2 CpG-DMR mCA level in D2-MSN(AM-seq)")
lm_FIt_WG_1 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_WG_1)
#lm_FIT_10_WG1=lm(mCA???poly(mCG ,10),data=dat_WG_1)


dat_WG_2 <- read.table("wgbs_D2_Rep2.CG.CA.txt")
colnames(dat_WG_2) <- c("EM1","mCG","mCA")
p6=ggplot(dat_WG_2,aes(x=mCG,y=mCA))+geom_point(alpha=0.05) + 
  geom_smooth(method='lm', formula= y~x+I(x^2)+I(x^3)+I(x^4),size=2,col="skyblue2")+xlab("D2 CpG-DMR meCG in D2-MSN(AM-seq)")+ylab("D2 CpG-DMR mCA level in D2-MSN(AM-seq)")
lm_FIt_WG_2 <- lm(mCA~mCG+I(mCG^2)+I(mCG^3)+I(mCG^4),data=dat_WG_2)
#lm_FIT_10_WG2=lm(mCA???poly(mCG ,10),data=dat_WG_2)

library(ggpubr)
theme_set(theme_pubr())
library("gridExtra")

#figure <- ggarrange(plotlist =list(p2,p3,p4),labels=c("A","B","C"),ncol=2,nrow=1)

# figure1 <- grid.arrange( p1.2, p2.2, p3.2,p4.2,p5.2,p6.2,layout_matrix = matrix(c(1, 2, 3, 4,5,6), nrow = 6))
figure1 <- grid.arrange( p1, p3,p5,layout_matrix = matrix(c(1,2,3), nrow =1))

ggsave("mCG_mcA.pdf",plot=figure1, device="pdf",width=12, height=4)
