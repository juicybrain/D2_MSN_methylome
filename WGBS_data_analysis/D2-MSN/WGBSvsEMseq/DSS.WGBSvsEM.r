#!/usr/bin/Rscript
library(DSS)
library(bsseq)





         dat_D2_WGBS1 = read.table("../DSS_GT_3/wgbs_D2_Rep1.DSS_GT3.DSS", header=F)
         colnames(dat_D2_WGBS1) <- c("chr",     "pos",     "N",       "X")

         dat_D2_WGBS2 = read.table("../DSS_GT_3/wgbs_D2_Rep2.DSS_GT3.DSS", header=F)
         colnames(dat_D2_WGBS2) <- c("chr",     "pos",     "N",       "X")

         dat_D2_EM1 = read.table("../DSS_GT_3/EM_D2_Rep1.DSS_GT3.DSS", header=F)
         colnames(dat_D2_EM1) <- c("chr",     "pos",     "N",       "X")

         dat_D2_EM2 = read.table("../DSS_GT_3/EM_D2_Rep2.DSS_GT3.DSS", header=F)
         colnames(dat_D2_EM2) <- c("chr",     "pos",     "N",       "X")




BSobj = makeBSseqData( list(dat_D2_WGBS1,dat_D2_WGBS2,dat_D2_EM1,dat_D2_EM2),c("WGBS1","WGBS2","EM1","EM2") )

#Treatment = factor(rep(c("R1","R2"),2))
#Pair = factor(rep(1:2,each=2))
#design = data.frame(Treatment, Pair)


#    DMLfit = DMLfit.multiFactor(BSobj, design, formula = ~ Treatment + Pair)
#    DMLfit$X
#    dmlTest.Treatment = DMLtest.multiFactor(DMLfit, term="Treatment")
#    dmlTest.Pair = DMLtest.multiFactor(DMLfit, term="Pair")
    # DMRs per Treatment
#           ix=sort(DMLtest.Treatment[,"pvals"], index.return=TRUE)$ix
#            DMR_Treatment=callDMR(DMLtest.Treatment, p.threshold=0.05)
#            write.table(DMR_Treatment,"R1vsR2_per_treatment.txt",sep="\t",quote=F,row.name=FALSE)
    # DMRs per pair
#            ix=sort(DMLtest.Pair[,"pvals"], index.return=TRUE)$ix
#            DMR_Pair=callDMR(DMLtest.Pair, p.threshold=0.05)
#            write.table(DMR_Pair,"R1vsR2_per_pair.txt",sep="\t",quote=F,row.name=FALSE)



dmlTest = DMLtest(BSobj, group1=c("WGBS1","WGBS2"), group2=c("EM1","EM2"),smoothing=F)
# dmlTest.sm = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"),
#                     smoothing=TRUE)
# call DML
#    dmls = callDML(dmlTest,delta=0.05, p.threshold=0.05)

    dmls = callDML(dmlTest, delta=0.1, p.threshold=0.001)
    write.table(dmls,"dml_WGBSvsEM_delta_0.001_nosmooth.txt",sep="\t",row.name=FALSE,quote=F)
#   system("cat dml_WGBSvsEM_delta_0.05_nosmooth.txt |sed '1d'|sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\t\"}{print $1,$2-1,$2,$5}'> dml_WGBSvsEM_delta_0.05_nosmooth.bed")


    dmrs = callDMR(dmlTest,delta=0.05,minlen=25,minCG=3)
#   dmrs = callDMR(dmlTest, delta=0.2, p.threshold=1e-5,minlen=250, minCG=4, dis.merge=250, pct.sig=0.5)
    write.table(dmrs,"dmr_WGBSvsEM_0.051e-5_nosmooth.txt",sep="\t",quote=F,row.name=FALSE)
    system("cat dmr_WGBSvsEM_0.051e-5_nosmooth.txt| sed '1d'| sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\t\"}{print $1,$2,$3}' > dmr_WGBSvsEM_0.051e-5_nosmooth.bed")
