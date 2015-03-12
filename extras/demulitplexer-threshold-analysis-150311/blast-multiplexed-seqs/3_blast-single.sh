# Constant
db=../blast-database/cell-line-db
outpref=3out_blast-output

# Derive
fastain=$1
idx=${fastain##*-}
idx=${idx%%.gz}
outname="$outpref"_chunck-"$idx".xml.gz

# Run
gzcat $fastain | blastn -db $db -outfmt 5 -num_alignments 2 | gzip -c > $outname
