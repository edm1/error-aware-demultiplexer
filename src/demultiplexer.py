#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
#

# import gzip
# import sys
# import os
# from shutil import rmtree
# import inspect
# from operator import itemgetter
# from src.bio_file_parsers import fastq_parser, phred_score_dict, write_fastq
# from Bio import SeqIO
# from src.probabilistic_seq_match import base_prob, sequences_match_prob
# import glob

import src.probabilistic_seq_match as seqprob
import sys
import os
from shutil import rmtree
import glob
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

def run(args):

    # Precompute base probabilities for phredscores up to 50
    base_prob_precompute = {}
    for score in range(1, 51):
        base_prob_precompute[score] = seqprob.base_prob(score)

    # Check that the multiplexed path exists
    multiplexed_dir = os.path.join(args.inDir, "multiplexed")
    if not os.path.exists(multiplexed_dir):
        sys.exit("Directory '<inDir>/multiplexed' does not exist. Re-run with"
                 " different <inDir>")

    # Create out directory
    out_dir = "demultiplexed"
    if args.uniqID != None:
        out_dir += "_{0}".format(args.uniqID)
    out_dir = os.path.join(args.inDir, out_dir)
    create_folder(out_dir)

    # Initiate multiplexed class
    multiplexed = Multiplex(multiplexed_dir)
    for out in multiplexed.iterate():
        rec = out[0][0]
        print(dir(rec))
        print(rec.letter_annotations['phred_quality'])
        break

    # Initiate sample sheet and read possible indexes
    sampleSheet = SampleSheet(args.sampleSheet)
    sampleSheet.parse(args.indexQual)
    for sample in sampleSheet.sample_indexes:
        rec = sampleSheet.sample_indexes[sample][0]
        print(rec.seq)
        print(rec.letter_annotations['phred_quality'])

    return 0

class SampleSheet:
    # Class to hold the sample sheet and retrieve indexes from it.

    def __init__(self, path):
        self.path = path

    def parse(self, index_qual):
        """ Parses the sample sheet to retrieve the indexes for each sample.
        """
        sample_indexes = {}
        with open(self.path, 'r') as in_h:
            # Skip to line after [Data]
            line = in_h.readline()
            while not line.startswith('[Data]'):
                line = in_h.readline()
            # Get header
            header = in_h.readline().rstrip().lower().split(',')
            col_ind = dict(zip(header, range(len(header))))
            # Save whether it is dual indexed
            if "index2" in col_ind.keys():
                self.is_dualindexed = True
            else:
                self.is_dualindexed = False
            # Get indexes
            for line in in_h:
                # Break if EOF
                if line.strip() == "":
                    break
                # Get info
                parts = line.rstrip().split(',')
                sample_name = parts[col_ind['sample_name']]
                # Get first index
                index1 = parts[col_ind['index']]
                sample_indexes[sample_name] = [index1]
                # Get second index
                if self.is_dualindexed:
                    index2 = parts[col_ind['index2']]
                sample_indexes[sample_name].append(index2)
        
        # Convert indexes to seqIO seqRecords
        self.sample_indexes = self.convert_index_to_seqRecord(sample_indexes,
            index_qual)

        return 0

    def convert_index_to_seqRecord(self, sample_indexes, index_qual):
        """ Converts each index sequence to a seqIO seqRecord.
        """
        # For each sample
        for sample in sample_indexes:
            # For each index
            for i in range(len(sample_indexes[sample])):
                raw_seq = sample_indexes[sample][i]
                # Convert to seqRecord
                sample_indexes[sample][i] = SeqRecord(
                    Seq(raw_seq, SingleLetterAlphabet()))
                # Add quality information
                sample_indexes[sample][i].letter_annotations['phred_quality'] \
                    = [index_qual] * len(raw_seq)
        return sample_indexes

class Multiplex:
    # Class for the folder of multiplexed reads + barcodes

    def __init__(self, folder):
        """ Make list of read and barcode files.
        """
        self.dir = folder
        # Get list of read and barcode paths
        self.read_paths = []
        self.barcode_paths = []
        for fastq in sorted(glob.glob(os.path.join(folder, "*.fastq*"))):
            if "barcode_" in os.path.split(fastq)[1]:
                self.barcode_paths.append(fastq)
            else:
                self.read_paths.append(fastq)
        # Save whether pairend
        if len(self.read_paths) == 1:
            self.is_pairend = False
        elif len(self.read_paths) == 2:
            self.is_pairend = True
        else:
            sys.exit("There must be 1 or 2 input read fastqs, not {0}".format(
                len(self.read_paths)))
        # Save whether dualindex
        if len(self.barcode_paths) == 1:
            self.is_dualindexed = False
        elif len(self.barcode_paths) == 2:
            self.is_dualindexed = True
        else:
            sys.exit("There must be 1 or 2 input barcode fastqs, not"
                     " {0}".format(len(self.barcode_paths)))
        return None

    def open_handles(self):
        """ Opens the file names for reading.
        """
        read_handles = [get_handle(filen, 'r') for filen in self.read_paths]
        barcode_handles = [get_handle(filen, 'r') for filen
                           in self.barcode_paths]
        return read_handles, barcode_handles

    def open_seqIO_iterators(self, read_handles, barcode_handles):
        """ Opens fastq iterators using biopythons SeqIO
        """
        # Open iterators for each handle
        read_iterators = [SeqIO.parse(handle, "fastq") for handle
                          in read_handles]
        barcode_iterators = [SeqIO.parse(handle, "fastq") for handle
                             in barcode_handles]
        return read_iterators, barcode_iterators

    def iterate(self):
        """ Loads the reads and barcode fastqs and yields 1 set at a time.
        """
        # Open handles
        read_handles, barcode_handles = self.open_handles()
        # Open iterators for each handle
        read_iterators, barcode_iterators = self.open_seqIO_iterators(
            read_handles, barcode_handles)
        print('here')

        # Iterate through records
        for r1_record in read_iterators[0]:
            
            # Get read records
            read_records = [r1_record]
            if self.is_pairend:
                read_records.append(next(read_iterators[1]))
            # Get barcode records
            barcode_records = [next(barcode_iterators[0])]
            if self.is_dualindexed:
                barcode_records.append(next(barcode_iterators[1]))

            # Check that they all have the same title
            titles = [record.id for record in read_records + barcode_records]
            if len(set(titles)) > 1:
                sys.exit('Reads and/or barcodes are not in sync\n'
                         '{0}'.format(titles))

            yield [read_records, barcode_records]

        # Close handles
        for handle in read_handles + barcode_handles:
            handle.close()


def create_folder(folder):
    """ Check out folder exists and create a new one.
    """
    # Check if it exists
    if not os.path.exists(folder):
        os.makedirs(folder)
    if os.path.exists(folder):
        response = input('{0} exists. Would you like to overwrite it? [y/n] '.format(folder))
        if response == 'y':
            rmtree(folder)
        else:
            sys.exit()
    os.makedirs(folder)
    return folder


def get_handle(filen, rw):
    """ Returns file handle using gzip if file ends in .gz
    """
    if filen.split('.')[-1] == 'gz':
        return gzip.open(filen, rw)
    else:
        return open(filen, rw)

"""
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
            print('Done: {0:,}'.format(c))
        c += 1

    print("Finished!")

    # Close all handles
    indexes.close_out_handles()
    multiplexed_reads.close_handles()
"""

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

        for sample_name in list(self.index_dict.keys()) + ['not_assigned']:

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

def phred_score_dict(offset):
    """ Creates a dict of phred score values
    """
    ascii_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    offset_string = ascii_string[offset - 33:]

    phred_dict = {}
    for i in range(len(offset_string)):
        phred_dict[offset_string[i]] = float(i)

    return phred_dict
