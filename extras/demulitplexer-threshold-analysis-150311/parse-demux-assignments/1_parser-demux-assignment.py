#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import sys
import glob
import gzip
import os

def main():

    # Options
    # demux_dir = "../data-from-jack/M01520_93_000000000-A72MP/demultiplexed-head/demultiplexed_prob0.05"
    demux_dir = sys.argv[1]

    # Get run name
    runname = os.path.split(demux_dir)[1]
    outname = "{0}_assingments.tsv.gz".format(runname)

    # Get list of cell line folders
    cellline_folders = glob.glob(os.path.join(demux_dir, "*"))

    with gzip.open(outname, "w") as out_h:
        for cellline_path in cellline_folders:
            cellline = os.path.split(cellline_path)[1]
            # Get name of fastq
            in_fastq = glob.glob(os.path.join(cellline_path,
                                              "*.R1.fastq.gz"))[0]
            # Parse fastq
            with gzip.open(in_fastq, "r") as in_h:
                for line in in_h:
                    if line.startswith("@M"):
                        readname = line.rstrip().lstrip("@")

                        # Output readname and cellline
                        lineout = [readname, cellline]
                        out_h.write("\t".join(lineout) + "\n")
                    

    return 0
   
if __name__ == '__main__':
    main()

 
