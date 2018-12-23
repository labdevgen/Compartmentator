java -jar ~/.local/bin/juicer_tools.1.9.8_jcuda.0.8.jar dump oe KR http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/Aste_VARYA_comb_allValidPairs.hic 2L 2L BP 25000 > Aste_VARYA_comb_allValidPairs.2L.oe
awk '($2-$1)<=10000000 {if (NR>1) print }' Aste_VARYA_comb_allValidPairs.2L.oe | head
awk '($2-$1)<=10000000 {if (NR>1) print }' Aste_VARYA_comb_allValidPairs.2L.oe > Aste_VARYA_comb_allValidPairs.2L.oe.10MB


awk '($2-$1)<=10000000 {if (NR>1) print }' Aste_VARYA_comb_allValidPairs.2L.10k.oe > Aste_VARYA_comb_allValidPairs.2L.10k.oe.10MB