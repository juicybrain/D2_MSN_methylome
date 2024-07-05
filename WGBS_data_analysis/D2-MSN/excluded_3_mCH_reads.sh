for spl in  D2_WGBS
do

    samtools view ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | awk 'BEGIN{OFS="\t"}{split($14,a,":")}{if(a[3] ~ /[a-z]/); else if(gsub(/[XH]/, "", a[3])>2)print $0}' | cut -f 1  > ${spl}\.new.SE.list

    samtools view -h ./${spl}\.R1_bismark_bt2_pe.deduplicated.bam | grep -vf ${spl}\.new.SE.list | samtools view -bS -o ${spl}\_filter.bam

#   samtools sort -n ${spl}\.bam > ${spl}\.ns.bam

#   bismark_methylation_extractor --CX --buffer_size 16G  --gzip --parallel 8 -p --bedGraph ${spl}\_filter.bam

done
