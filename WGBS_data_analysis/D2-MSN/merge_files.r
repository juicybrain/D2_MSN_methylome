#!/usr/bin/Rscript



dat_all = read.table(paste0("D2M_WGBS_1",".fdr_fix.CC.DSS"),header=F)
colnames(dat_all) = c("chr","start", paste0("D1M_WGBS_1_CG", "me"),  paste0("D1M_WGBS_1", "N"))


for (spl in c("D1M_WGBS_2",  "D1M_EM_1", "D1M_EM_2" )){
    dat_tmp= read.table( paste0(spl, ".fdr_fix.CC.DSS"), header=F)
    colnames(dat_tmp) = c("chr", "start", paste0(spl, "me"), paste0(spl, "N"))
    dat_all = merge(dat_all, dat_tmp, by=c("chr","start"), all=F)
}


write.table(dat_all, "all_sample.CC.DSS", quote=F, sep="\t", row.names=F, col.names=T)

