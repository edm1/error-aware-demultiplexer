#fastqin=../data-from-jack/M01520_93_000000000-A72MP/raw_extract.subset.fastq.gz
fastqin=../data-from-jack/M01520_93_000000000-A72MP/raw_extract.1.fastq.gz
fastaout=1out_raw-reads.fasta.gz

echo "Converting..."
gzcat $fastqin | fastq_to_fasta -Q33 -n | gzip -c > $fastaout
echo "Done!"
