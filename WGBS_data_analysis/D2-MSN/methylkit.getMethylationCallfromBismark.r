#!/home/yli/bin/Rscript3.6

library("methylKit")

spl=c("","")

spl_name=c("EM_10_S22.pe.s.bam","EM_11_S23.pe.s.bam","EM_12_S24.pe.s.bam","EM_1_S13.pe.s.bam","EM_2_S14.pe.s.bam","EM_3_S15.pe.s.bam","EM_4_S16.pe.s.bam","EM_5_S17.pe.s.bam","EM_6_S18.pe.s.bam","EM_7_S19.pe.s.bam","EM_8_S20.pe.s.bam","EM_9_S21.pe.s.bam")

for (i in 1:12){
my.methRaw=processBismarkAln( location =spl[i], mincov=3,sample.id=spl_name[i], assembly="mm10",read.context="CpG", save.folder=getwd())
}
