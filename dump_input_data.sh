hic="http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/Aste_VARYA_comb_allValidPairs.hic";
juicer="$HOME/.local/bin/juicer_tools.1.9.8_jcuda.0.8.jar";
res="10000";

out=$(echo $hic | rev | cut -d "/" -f1 | rev)".$res.oe"; \
rm $out;
for chr in {"2L","2R","3L","3R","X"}; do
    echo $chr;
    java -jar $juicer dump oe KR $hic $chr $chr BP $res | \
    awk -v chr="$chr" 'NR>1 {OFS="\t"; print chr,$0}' >> $out;
done;
awk '($3-$2)<=10000000' $out > $out.10MB