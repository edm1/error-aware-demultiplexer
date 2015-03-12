fastapattern=2out_raw-reads_chunck-*.gz

parallel -j 12 bash 3_blast-single.sh ::: $fastapattern
