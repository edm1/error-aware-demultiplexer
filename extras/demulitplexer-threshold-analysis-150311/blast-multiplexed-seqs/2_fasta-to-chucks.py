#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import subprocess
import gzip

def main():

    in_fasta = "1out_raw-reads.fasta.gz"
    out_pref = "2out_raw-reads_chunck-"
    num_chuncks = 4

    # Count seqs
    total_seqs = count_num_seqs(in_fasta)
    per_chuck = int(total_seqs / num_chuncks) + 1
    lines_per_file = per_chuck * 2

    # Use split command
    cmd = "gzcat {1} | split -l {0} - {2}".format(lines_per_file, in_fasta, out_pref)
    subprocess.call(cmd, shell=True)

    # Gzip all chuncks
    cmd = "gzip {0}*".format(out_pref)
    subprocess.call(cmd, shell=True)

    return 0

def count_num_seqs(filen):
    """ Counts the number of sequences in a fasta
    """
    num = 0
    with gzip.open(filen, "r") as in_h:
        for line in in_h:
            if line.startswith(">"):
                num += 1
    return num


if __name__ == '__main__':
    main()

 
