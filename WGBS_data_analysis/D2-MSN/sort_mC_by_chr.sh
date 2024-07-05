    for i in {1..19}
    do
    cat F_N_meCA.bedgraph  F_N_meCA.part2.bedgraph  | grep chr${i}$'\t' |sort -k1,1 -k2,2n > split_and_sort/F_N_chr${i}\.s.bedgraph &
    cat F_P_meCA.bedgraph  F_P_meCA.part2.bedgraph  | grep chr${i}$'\t' |sort -k1,1 -k2,2n > split_and_sort/F_P_chr${i}\.s.bedgraph &
    cat M_N_meCA.bedgraph  M_N_meCA.part2.bedgraph  | grep chr${i}$'\t' |sort -k1,1 -k2,2n > split_and_sort/M_N_chr${i}\.s.bedgraph &
    cat M_P_meCA.bedgraph  M_P_meCA.part2.bedgraph  | grep chr${i}$'\t' |sort -k1,1 -k2,2n > split_and_sort/M_P_chr${i}\.s.bedgraph

    done
