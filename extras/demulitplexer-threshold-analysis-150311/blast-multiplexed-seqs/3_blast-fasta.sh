fastain=2out_raw-reads_chunck-*
db=../blast-database/cell-line-db
outpref=3out_blast-output

for filen in $fastain; do
	idx=${filen##*-}
	idx=${idx%%.gz}
	echo $idx
	outname="$outpref"_chunck-"$idx".xml.gz
	echo "gzcat $filen | blastn -db $db -outfmt 5 -num_alignments 2 | gzip -c > $outname "
done
