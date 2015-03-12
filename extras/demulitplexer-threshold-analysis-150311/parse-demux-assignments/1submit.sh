indir=/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/data-from-jack/M01520_93_000000000-A72MP
script=/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/parse-demux-assignments/1_parser-demux-assignment.py

for d in $(ls -d $indir/demult*); do
    python $script $d
done