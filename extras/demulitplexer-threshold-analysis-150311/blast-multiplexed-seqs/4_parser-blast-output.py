#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from Bio.Blast import NCBIXML
import sys
import glob
import gzip

def main():

    in_xml_pattern = "3out_blast-output_chunck-*.xml.gz"
    out_name = "4out_blast-output-parsed.tsv"

    # Open out handle
    with open(out_name, "w") as out_h:

        # Write header
        line = ["query", "best_hit", "best_score", "next_best_score"]
        out_h.write("\t".join(line) + "\n")

        # Parse blast xmls
        for in_xml in glob.glob(in_xml_pattern):
            with gzip.open(in_xml, "r") as xml:
                
                # Create parser
                blast_records = NCBIXML.parse(xml)

                # Process queries
                for query in blast_records:
                        
                    # Get results and write to output
                    res = find_best_hit(query)
                    line = [str(x) for x in res]
                    out_h.write("\t".join(line) + "\n")

    return 0

def find_best_hit(query):
    """ Takes a blast query record and returns the hit and bit score.
        Also returns bit score of next best hit if any.
    """

    query_id = query.query
    best_hit = "No hit"
    best_score = None
    next_score = None

    # If there are no alignments then skip
    if not len(query.alignments) == 0:
        
        # Get hit of best alignment
        best_hit = query.alignments[0].title
        best_score = query.alignments[0].hsps[0].bits
        
        # If there is a 2nd hit then get bit score of that
        if len(query.alignments) > 1:
            next_score = query.alignments[1].hsps[0].bits

    return query_id, best_hit, best_score, next_score

if __name__ == '__main__':
    main()

 
