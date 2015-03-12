fastqin=../data-from-jack/M01520_93_000000000-A72MP/raw_extract.subset.fastq.gz
#fastqin=../data-from-jack/M01520_93_000000000-A72MP/raw_extract.1.fastq.gz
fastaout=1out_raw-reads.fasta.gz

gzcat $fastqin | fastq_to_fasta -n | gzip -c > $fastaout
