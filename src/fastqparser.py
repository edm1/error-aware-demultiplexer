#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Contains simple parsers for fasta and fastq files, taken directly from the
# biopython source code.
#

import string

def wrap(string, length):
    """ Yield successive length-sized chunks from string.
    """
    for i in xrange(0, len(string), length):
        yield string[i:i + length]

def phred_score_dict(offset):
    """ Creates a dict of phred score values
    """
    phred_dict = {}
    for letter in string.printable:
        phred_dict[letter] = float(ord(letter) - offset)
    return phred_dict

def fastqIterator(handle):
    """Iterate over Fastq records 
    """
    # We need to call handle.readline() at least four times per record,
    # so we'll save a property look up each time:
    handle_readline = handle.readline

    # Skip any text before the first record (e.g. blank lines, comments?)
    while True:
        line = handle_readline()
        if not line:
            return  # Premature end of file, or just empty?
        if line[0] == "@":
            break
        if isinstance(line[0], int):
            raise ValueError("Is this handle in binary mode not text mode?")

    while line:
        if line[0] != "@":
            raise ValueError(
                "Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()
        # Will now be at least one line of quality data - in most FASTQ files
        # just one line! We therefore use string concatenation (if needed)
        # rather using than the "".join(...) trick just in case it is multiline:
        seq_string = handle_readline().rstrip()
        # There may now be more sequence lines, or the "+" quality marker line:
        while True:
            line = handle_readline()
            if not line:
                raise ValueError("End of file without quality information.")
            if line[0] == "+":
                # The title here is optional, but if present must match!
                second_title = line[1:].rstrip()
                if second_title and second_title != title_line:
                    raise ValueError("Sequence and quality captions differ.")
                break
            seq_string += line.rstrip()  # removes trailing newlines
        # This is going to slow things down a little, but assuming
        # this isn't allowed we should try and catch it here:
        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)

        # Will now be at least one line of quality data...
        quality_string = handle_readline().rstrip()
        # There may now be more quality data, or another sequence, or EOF
        while True:
            line = handle_readline()
            if not line:
                break  # end of file
            if line[0] == "@":
                # This COULD be the start of a new sequence. However, it MAY just
                # be a line of quality data which starts with a "@" character.  We
                # should be able to check this by looking at the sequence length
                # and the amount of quality data found so far.
                if len(quality_string) >= seq_len:
                    # We expect it to be equal if this is the start of a new record.
                    # If the quality data is longer, we'll raise an error below.
                    break
                # Continue - its just some (more) quality data.
            quality_string += line.rstrip()

        if seq_len != len(quality_string):
            raise ValueError("Lengths of sequence and quality values differs "
                             " for %s (%i and %i)."
                             % (title_line, seq_len, len(quality_string)))

        # Convert into a fastq record
        record = Fastq(title_line, seq_string, quality_string)

        # Return the record and then continue...
        yield record
    raise StopIteration

def fastqWriter(record, handle):
    """ Simple fastq writer.
    """
    out = []
    # Add @ to header
    out.append("@{0}".format(record.id))
    # Add sequence
    out.append(record.seq)
    # Add +
    out.append("+")
    # Add qual
    out.append(record.qual_str)
    # Write to handle
    out_line = "\n".join(out) + "\n"
    handle.write(out_line.encode("utf-8"))
    return 0


class Fastq:
    # Class to hold a single fastq record

    def __init__(self, title_line, seq_string, quality_string):
        self.id = title_line
        self.seq = seq_string
        self.qual_str = quality_string
        self.qual_prob = None

    def qual_to_prob(self, base_prob_precompute):
        """ Converts the quality string to a list of base probabilities.
        """
        self.qual_prob = [base_prob_precompute[x] for x in self.qual_str]
        return 0



