#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
#

import argparse
import gzip
import sys
import os
from shutil import rmtree
import inspect
from operator import itemgetter
from src.bio_file_parsers import fastq_parser, phred_score_dict, write_fastq
from src.probabilistic_seq_match import base_prob, sequences_match_prob
import glob

class Indexes:
    # Class to parse and hold the indexes.

    def __init__(self):
        self.index_dict = None
        self.dual_end = False
        self.pair_end = False
        self.demux_handles = None
        self.barcode_handles = None

    def open_out_handles(self, results_dir):
        """ For each sample name open an write handle.
        """
        # For each sample name create a output file handle
        self.demux_handles = {}
        self.barcode_handles = {}

        for sample_name in self.index_dict.keys() + ['not_assigned']:

            # Open handles for the read 1
            reads_name = os.path.join(results_dir, '{0}.R1.fastq.gz'.format(sample_name.replace(' ', '_')))
            self.demux_handles[sample_name] = [get_handle(reads_name, 'w')]
            # If pair-end then append another handle
            if self.pair_end == True:
                reads_name = os.path.join(results_dir, '{0}.R2.fastq.gz'.format(sample_name.replace(' ', '_')))
                self.demux_handles[sample_name].append(get_handle(reads_name, 'w'))

            # Open for barcode 1
            barcode_name = os.path.join(results_dir, '{0}.barcode_1.fastq.gz'.format(sample_name.replace(' ', '_')))
            self.barcode_handles[sample_name] = [get_handle(barcode_name, 'w')]
            # If dual indexed then append
            if self.dual_end == True:
                barcode_name = os.path.join(results_dir, '{0}.barcode_2.fastq.gz'.format(sample_name.replace(' ', '_')))
                self.barcode_handles[sample_name].append(get_handle(barcode_name, 'w'))

        return 0

    def num_indexes(self):
        if self.dual_end == True:
            return 2
        else:
            return 1

    def close_out_handles(self):
        """ Closes all of the open out handles.
        """
        for handle_dict in [self.demux_handles, self.barcode_handles]:
            for sample_name in handle_dict:
                for handle in handle_dict[sample_name]:
                    handle.close()
        return 0

    def load_samplesheet_indexes(self, filen, index_ascii):
        """ Will find check if one of two indexes are used then load them.
        """
        with open(filen, 'r') as in_h:
            # Skip to line after [Data]
            line = in_h.readline()
            while not line.startswith('[Data]'):
                line = in_h.readline()
            # Get header
            header = in_h.readline().rstrip().lower().split(',')
            col_ind = dict(zip(header, range(len(header))))
            # Get indexes
            self.index_dict = {}
            for line in in_h:
                # Break if EOF
                if line.strip() == "":
                    break
                # Get info
                parts = line.rstrip().split(',')
                sample_name = parts[col_ind['sample_name']]
                # Get first index
                index1 = parts[col_ind['index']]
                self.index_dict[sample_name] = [(index1, index_ascii * len(index1))]
                # Get second index
                if "index2" in col_ind.keys():
                    index2 = parts[col_ind['index2']]
                    self.index_dict[sample_name].append((index2, index_ascii * len(index1)))
        # Save whether it is dual end
        if "index2" in col_ind.keys():
            self.dual_end = True

class Record:
    # Class to hold a single fastq record

    def __init__(self, record):
        """ record is a tuple/list in form (title, seq, qual)
        """
        self.title = record[0]
        self.seq = record[1]
        self.qual_string = record[2]

class Reads:
    # Class to iterate over the read/barcodes and return a set of each

    def __init__(self, in_dir):
        """ Searches the input directory for reads and barcode file names.
        """
        self.read_filenames = []
        self.barcode_filenames = []

        # Find read file names
        c = 1
        while True:
            read_files = glob.glob(os.path.join(in_dir, "*.{0}.fastq.gz".format(c)))
            if len(read_files) > 0:
                self.read_filenames.extend(read_files)
                c += 1
            else:
                break
        # Find the barcode file names
        c = 1
        while True:
            barcode_files = glob.glob(os.path.join(in_dir, "*.barcode_{0}.fastq.gz".format(c)))
            if len(barcode_files) > 0:
                self.barcode_filenames.extend(barcode_files)
                c += 1
            else:
                break

    def open_handles(self):
        """ Opens the file names for reading.
        """
        self.read_handles = [get_handle(filen, 'r') for filen in self.read_filenames]
        self.barcode_handles = [get_handle(filen, 'r') for filen in self.barcode_filenames]

    def close_handles(self):
        """ Closes open handes
        """
        for handle_list in [self.read_handles, self.barcode_handles]:
            for handle in handle_list:
                handle.close()

    def iterate(self):
        """ Loads the reads and barcode fastqs and yields 1 set at a time.
        """
        # Open iterators for each handle
        read_iterators = [fastq_parser(handle) for handle in self.read_handles]
        barcode_iterators = [fastq_parser(handle) for handle in self.barcode_handles]

        # Iterate trhough records in 1st read file
        for r1_record in read_iterators[0]:

            # Get read records
            read_records = []
            # Append 1st read record
            read_records.append(Record(r1_record))
            # Append subsequent
            for iterator in read_iterators[1:]:
                read_records.append(Record(iterator.next()))

            # Get barcode records
            barcode_records = [Record(iterator.next()) for iterator in barcode_iterators]

            # Check that they all have the same title
            titles = [record.title.split(" ")[0] for record in read_records + barcode_records]
            if len(set(titles)) > 1:
                sys.exit('Reads and/or barcodes are not in sync\n{0}'.format(titles))

            yield [read_records, barcode_records]

def main():

    # Get root of project directory
    global root_dir
    root_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    # Parse args
    args = parse_arguments()

    # Create a dictionary of phred scores
    phred_dict = phred_score_dict(args.PhredOffset)
    phred_dict_inv = {int(v):k for k, v in phred_dict.items()}

    # Precompute base probabilities
    base_prob_precompute = {}
    for key in phred_dict:
        base_prob_precompute[key] = base_prob(phred_dict[key])

    # Create out directory
    out_dir = create_results_folder(args.InputTempDir, args.ResultsDir, args.OutID)

    # Load reads/barcodes
    multiplexed_reads = Reads(args.InputTempDir)
    # Open file handles
    multiplexed_reads.open_handles()

    # Load indexes from sample sheet
    indexes = Indexes()
    indexes.load_samplesheet_indexes(args.SampleSheet, phred_dict_inv[args.IndexQual])
    # Update indexes class with whether pe or se
    if len(multiplexed_reads.read_filenames) == 2:
        indexes.pair_end = True
    # Create output handles in results folder for each sample name
    indexes.open_out_handles(out_dir)

    # Check number of barcode files matches number of indexes
    if not len(multiplexed_reads.barcode_filenames) == indexes.num_indexes():
        sys.exit('Incorrect number of indexes in sample sheet. Exiting.')

    # Iterate over each barcode/read set and look for matches in indexes
    c = 1
    for read_records, barcode_records in multiplexed_reads.iterate():

        # Find best index match
        b1_header, sample, prob = match_barcode_to_indexes(barcode_records, indexes,
                                                           base_prob_precompute, args.MinProb)
        if sample == None:
            sample = 'not_assigned'

        # Append probability to the (b1) header
        b1_header = "{0} {1}".format(b1_header, prob)

        # Write fastqs
        read_records[0]
        for i in range(len(read_records)):
            write_fastq(indexes.demux_handles[sample][i],
                        read_records[i].title,
                        read_records[i].seq,
                        read_records[i].qual_string)

        # Write barcode fastqs
        for i in range(len(barcode_records)):
            write_fastq(indexes.barcode_handles[sample][i],
                        barcode_records[i].title,
                        barcode_records[i].seq,
                        barcode_records[i].qual_string)

        # Update progress
        if c % 10000 == 0:
            print 'Done: {0:,}'.format(c)
        c += 1

    print "Finished!"

    # Close all handles
    indexes.close_out_handles()
    multiplexed_reads.close_handles()


def create_results_folder(input_temp_dir, results_dir, out_id):
    """ Check out folder exists and create a new one.
    """

    out_dir = os.path.join(results_dir, input_temp_dir.rstrip('/').split('/')[-1])
    out_dir += '_{0}'.format(out_id)
    # Check if it exists
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    if os.path.exists(out_dir):
        response = raw_input('{0} exists. Would you like to overwrite it? [y/n] '.format(out_dir))
        if response == 'y':
            rmtree(out_dir)
        else:
            sys.exit()
    os.makedirs(out_dir)
    return out_dir


def match_barcode_to_indexes(barcode_records, indexes, base_prob_precompute, min_prob):
    """ For the barcode pair, caluclates probability of a match against each set
        of indexes
    """

    index_probs = {}

    for sample_name in indexes.index_dict:

        index_records = indexes.index_dict[sample_name]

        # Calculate the match probability for barcode 1
        b1_prob = sequences_match_prob(index_records[0][0], index_records[0][1],
                                       barcode_records[0].seq, barcode_records[0].qual_string,
                                       base_prob_precompute, 0)

        # If dual indexed calc prob of other match
        b2_prob = None
        if indexes.dual_end == True:
            # Skip if already below the threshold
            if b1_prob >= min_prob:
                b2_prob = sequences_match_prob(index_records[1][0], index_records[1][1],
                                               barcode_records[1].seq, barcode_records[1].qual_string,
                                               base_prob_precompute, 0)
            else:
                b2_prob = b1_prob

        # Caluclate combined probability
        if b2_prob != None:
            overall_prob = b1_prob * b2_prob
        else:
            overall_prob = b1_prob
        index_probs[sample_name] = overall_prob

    sorted_probs = sorted(index_probs.iteritems(), key=itemgetter(1), reverse=True)

    # Return header, sample, prob
    header = barcode_records[0].title
    if sorted_probs[0][1] > min_prob:
        return header, sorted_probs[0][0], sorted_probs[0][1]
    else:
        return header, None, sorted_probs[0][1]


def get_handle(filen, rw):
    """ Returns file handle using gzip if file ends in .gz
    """
    if filen.split('.')[-1] == 'gz':
        return gzip.open(filen, rw)
    else:
        return open(filen, rw)


def parse_arguments():
    """ Load the arguments required to run IlluminaBasecallsToFastq.
    """

    parser = argparse.ArgumentParser(description='Run IlluminaBasecallsToFastq using picard.')

    # Required args
    parser.add_argument('--InputTempDir',
                        metavar='<dir>',
                        type=str,
                        required=True,
                        help='Directory containing files output by 1_run_IlluminaBasecallsToFastq.py')
    parser.add_argument('--SampleSheet',
                        metavar='<SampleSheet.csv>',
                        type=str,
                        required=True,
                        help='MiSeq output SampleSheet.csv file')
    parser.add_argument('--OutID',
                        metavar='<str>',
                        type=str,
                        required=True,
                        help='MiSeq output SampleSheet.csv file')

    # Optional arguments
    parser.add_argument('--ResultsDir',
                        metavar='<dir>',
                        type=str,
                        required=False,
                        default=os.path.join(root_dir, 'results'),
                        help='Results directory. (./results)')
    parser.add_argument('--MinProb',
                        metavar='<float>',
                        type=float,
                        required=False,
                        default=0.05,
                        help='Minimum probability of a match else discard. (0.05)')
    parser.add_argument('--PhredOffset',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=33,
                        help='FASTQ phred score offset (33)')
    parser.add_argument('--IndexQual',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=30,
                        help='Phred-score given to barcode indexes (30)')
    parser.add_argument('--version', action='version', version='v0.2')


    return parser.parse_args()

if __name__ == '__main__':
    main()
