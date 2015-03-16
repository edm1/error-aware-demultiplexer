#!/usr/bin/env bash
script=2_process-single.R
filein_pattern=1out_blast_demux_merge_probPROB.tsv.gz

for prob in 0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 0.7 0.9; do
	inname=${filein_pattern/PROB/$prob}
	outname=2out_prob"$prob"_table.txt
	cmd="Rscript "$script" "$inname" "$outname""
	echo $cmd
done | parallel -j 16