for spl in `ls *fdr_fixed.CA.DSS.bed`
do

bedtools map -a /home/yli/ref/bed/windows/D1_WGBS_covered_regions.100k.s.f.bed  -b <(zcat ${spl} ) -c 7,8  -o sum > 100k_whole_genome/${spl}\.100k.bed &
done

