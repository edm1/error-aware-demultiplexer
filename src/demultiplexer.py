# -*- coding: utf-8 -*-
#
#
#

from src.probabilistic_seq_match import sequences_match_prob
from src.probabilistic_seq_match import base_prob
import sys
import os
from shutil import rmtree
import glob
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from operator import itemgetter
import concurrent.futures as cf

def run(args):

    # Precompute base probabilities for phredscores up to 50
    base_prob_precompute = {}
    for score in range(1, 51):
        base_prob_precompute[score] = base_prob(score)

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

    # Initiate sample sheet and read possible indexes
    sampleSheet = SampleSheet(args.sampleSheet)
    sampleSheet.parse(args.indexQual)

    # Check that there are the same number of indexes in sample sheet and
    # multiplexed fastqs
    if sampleSheet.is_dualindexed != multiplexed.is_dualindexed:
        sys.exit("Error: Different number of indexes in sampleSheet and "
                 "multiplexed reads. Exiting!")

    # Open output class for each sample, and a not_assigned group
    sample_out = {}
    for sample in list(sampleSheet.sample_indexes.keys()) + ['not_assigned']:
        sample_out[sample] = Sample(sample, out_dir, multiplexed.is_pairend,
                                    multiplexed.is_dualindexed)

    # Version that doesnt use concurrent.futures
    c = 1
    for variables in futures_iterate_reads(multiplexed, sampleSheet,
                                           base_prob_precompute, args.minProb):
        output = futures_barcode_to_indexes(variables)
        # Unpack output
        ((read_records, barcode_records), sample, prob, _) = output
        # Write record to correct sample file
        sample_out[sample].write(read_records, barcode_records)
        # Update progress
        if c % 1000 == 0:
            print(c)
        c += 1


    """
    # Send each read/barcode record to futures to match up to sample
    with cf.ProcessPoolExecutor(max_workers=args.numCPU) as executor:
        c = 1
        # Map read/barcode records
        for output in executor.map(futures_barcode_to_indexes,
            futures_iterate_reads(multiplexed, sampleSheet,
                base_prob_precompute, args.minProb)):
            # Unpack output
            ((read_records, barcode_records), sample, prob, _) = output
            # Write record to correct sample file
            sample_out[sample].write(read_records, barcode_records)

            # Update progress
            if c % 1000 == 0:
                print(c)
            c += 1
    """

    return 0

def futures_iterate_reads(multiplexed, sampleSheet, base_prob_precompute,
                          min_prob):
    """ Returns an iterator that contains everything needed for futures.
    """
    for combined_record in multiplexed.iterate():
        yield (combined_record, sampleSheet, base_prob_precompute, min_prob)

def futures_barcode_to_indexes(variables):
    """ Compares the reads barcodes to sample indexes and returns matching
        sample name.
    """
    # Unpack variables
    (combined_record, sampleSheet,
     base_prob_precompute, min_prob) = variables
    # Get barcode records
    _, barcode_records = combined_record
    # Find sample
    b1_header, sample, prob = match_barcode_to_indexes(barcode_records,
        sampleSheet, base_prob_precompute, min_prob)
    if sample == None:
        sample = 'not_assigned'
    # Append probability to barcode1 header
    b1_header = "{0}_{1}".format(b1_header, prob)
    # Change header
    combined_record[1][0].id = b1_header

    return combined_record, sample, prob, b1_header


def match_barcode_to_indexes(barcode_records, sampleSheet, base_prob_precompute,
                             min_prob):
    """ For the barcode pair, caluclates probability of a match against each set
        of indexes
    """

    index_probs = {}

    for sample_name in sampleSheet.sample_indexes:

        index_records = sampleSheet.sample_indexes[sample_name]

        # Calculate the match probability for barcode 1
        b1_prob = sequences_match_prob(index_records[0].seq,
            index_records[0].letter_annotations['phred_quality'],
            barcode_records[0].seq,
            barcode_records[0].letter_annotations['phred_quality'],
            base_prob_precompute, 0)

        # Do for second barcode if present
        if sampleSheet.is_dualindexed:
            # Skip if already below the threshold, else assign same prob as b1
            if b1_prob >= min_prob:
                b2_prob = sequences_match_prob(index_records[1].seq,
                           index_records[1].letter_annotations['phred_quality'],
                           barcode_records[1].seq,
                           barcode_records[1].letter_annotations['phred_quality'],
                           base_prob_precompute, 0)
            else:
                b2_prob = b1_prob

        # Caluclate combined probability
        if sampleSheet.is_dualindexed:
            overall_prob = b1_prob * b2_prob
        else:
            overall_prob = b1_prob

        index_probs[sample_name] = overall_prob

    sorted_probs = sorted(index_probs.items(), key=itemgetter(1),
                          reverse=True)

    # Return header, sample, prob
    header = barcode_records[0].id
    if sorted_probs[0][1] > min_prob:
        return header, sorted_probs[0][0], sorted_probs[0][1]
    else:
        return header, None, sorted_probs[0][1]


class Sample:
    # Class for each possible sample. 1) Holds the output directory for that
    # sample. 2) Opens handles. 3) Writes record to sample.
    
    def __init__(self, name, out_dir, is_pe, id_dual):
        self.read_paths = []
        self.barcode_paths = []
        self.read_handles = None
        self.barcode_handles = None
        
        # Create directory for sample
        name = name.replace(' ', '_')
        self.sample_dir = os.path.join(out_dir, name)
        create_folder(self.sample_dir)

        # Create read paths
        self.read_paths.append(os.path.join(self.sample_dir,
            '{0}.R1.fastq.gz'.format(name)))
        if is_pe:
            self.read_paths.append(os.path.join(self.sample_dir,
                '{0}.R2.fastq.gz'.format(name)))
        # Create barcode paths
        self.barcode_paths.append(os.path.join(self.sample_dir,
            '{0}.barcode_1.fastq.gz'.format(name)))
        if id_dual:
            self.barcode_paths.append(os.path.join(self.sample_dir,
                '{0}.barcode_2.fastq.gz'.format(name)))

    def open_handles(self):
        """ For the reads and barcodes, opens output handles.
        """
        self.read_handles = [get_handle(read_path, 'wt') for read_path
                             in self.read_paths]
        self.barcode_handles = [get_handle(barcode_path, 'wt') for barcode_path
                                in self.barcode_paths]
        return 0

    def write(self, read_records, barcode_records):
        """ Writes the demultiplexed read and barcode records to sample file.
        """
        # Open handles if not open
        if self.read_handles == None:
            self.open_handles()
        # Write read records
        for i in range(len(read_records)):
            SeqIO.write(read_records[i], self.read_handles[i], "fastq")
        # Write barcode records
        for i in range(len(barcode_records)):
            SeqIO.write(barcode_records[i], self.barcode_handles[i], "fastq")
        return 0

    def close_handles():
        """ Closes any open handles.
        """
        if self.read_handles != None:
            for handle in self.read_handles + self.barcode_handles:
                handle.close()
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
