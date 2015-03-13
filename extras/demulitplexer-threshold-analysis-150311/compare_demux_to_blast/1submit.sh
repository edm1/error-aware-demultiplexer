demux_files=/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/parse-demux-assignments/demultiplexed_prob*_assingments.tsv.gz
blast_file=/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/blast-multiplexed-seqs/4out_blast-output-parsed.tsv.gz

parallel -j 16 python 1_merge-blast-demux.py $blast_file {} ::: $demux_files
