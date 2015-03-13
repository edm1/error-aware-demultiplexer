#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import gzip
import pandas as pd
import re
import os

def main():

    # Options
    # in_blast = "../blast-multiplexed-seqs/4out_blast-output-parsed_head.tsv.gz"
    # in_demux = "../parse-demux-assignments/demultiplexed_prob0.05_assingments_head.tsv.gz"
    in_blast = sys.argv[1]
    in_demux = sys.argv[2]

    # Get prob name
    probname = re.search(r"_(prob.+?)_", os.path.split(in_demux)[1]).group(1)
    
    # Load dfs
    blast_df = pd.read_csv(in_blast, sep="\t", compression="gzip", header=0)
    demux_df = pd.read_csv(in_demux, sep="\t", compression="gzip", header=None)
    demux_df.columns = ["query", "demux_cell_line"]

    # Merge on query
    merge_df = pd.merge(blast_df, demux_df, how="left", on="query")

    # Output
    outname = "1out_blast_demux_merge_{0}.tsv.gz".format(probname)
    merge_df.to_csv(outname, sep="\t", compression="gzip", header=True,
                    index=False, na_rep="NA")


    print("Done!")

    return 0
   
if __name__ == '__main__':
    main()

 
