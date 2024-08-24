## sampling bed file from bam files

## shuffle bed file in mapable regions


for spl in `ls *bed`

do

seqtk subseq   /home/yli/ref/genome/male.mm10/male.mm10.fasta ${spl} > ${spl}\_peaks.fa


perl ~/workspace/script/dinucleotide_counts.pl  ${spl}\_peaks.fa > ${spl}\.report

done

