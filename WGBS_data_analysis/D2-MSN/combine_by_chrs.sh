
for spl in EM_D2_Rep1 EM_D2_Rep2  wgbs_D2_Rep2

do

cd $spl

for chr in chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrM chrX chrY
    do
cat ${spl}\_${chr}.CG.cov| sort -k1,1 -k2,2n >> ${spl}\.CpG.cov


cat ${spl}\_${chr}.CA.cov | sort -k1,1 -k2,2n >> ${spl}\.CpA.cov

cat ${spl}\_${chr}.CT.cov | sort -k1,1 -k2,2n >> ${spl}\.CpT.cov

cat ${spl}\_${chr}.CC.cov | sort -k1,1 -k2,2n >> ${spl}\.CpC.cov
    done
wait
gzip ${spl}\.CpG.cov &
gzip ${spl}\.CpA.cov &
gzip ${spl}\.CpT.cov &
gzip ${spl}\.CpC.cov &
cd ..

done