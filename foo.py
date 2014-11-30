#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Edward Mountjoy
# 
# Parses arguments and executes
#

import argparse
import os
import inspect
from src import basecalls2fastq
from src import demultiplexer
import sys

def main():
    """ Parses args and executes required script.
    """
    # Find scripts root directory
    root_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    # Parse the command line args
    args = parse_arguments(root_dir)
    args.rootDir = root_dir

    # Run required function
    args.func(args)

    print("\nFinished.")

    return 0

def parse_arguments(root_dir):
    """ Will parse the command line arguments arnd return the arg object.
    """

    # Create top level parser. TODO: add description, usage, etc
    parser = argparse.ArgumentParser('main')
    subparsers = parser.add_subparsers()

    # Create parser for the bcl2fastq (extracting reads from illumina folder)
    parser_b2f = subparsers.add_parser('bcl2fastq')
    # Required positional arguments
    parser_b2f.add_argument('baseCallDir', metavar='<baseCallDir>', type=str,
        help='Directory containing base call intensitites')
    parser_b2f.add_argument('runParamXML', metavar='<runParameters.xml>',
        type=str, help='runParameters.xml file')
    parser_b2f.add_argument('lane', metavar='<lane>',
        type=int, help='Lane number')
    # Optional arguments
    parser_b2f.add_argument('--outDir', '-o', metavar='<str>', type=str,
        default=os.path.join(root_dir, "output"),
        help='Location to create output files. (output)')
    parser_b2f.add_argument('--numCPU', '-p', metavar='<int>', type=int,
        default=1, help='Number of CPUs to use. (1)')
    parser_b2f.add_argument('--readsPerTile', '-r', metavar='<int>', type=int,
        default=120000, help=('Max number of reads in RAM per tile, reduce if '
                              'you have problems with memory. (120000)'))
    parser_b2f.add_argument('--MaxInRam', '-m', metavar='<int>', type=int,
        default=500000, help=('Maximum number of records that are stored in the'
                              ' RAM. (500000)'))
    parser_b2f.add_argument('--JavaRAM', '-mem', metavar='<int>', type=int,
        default=2, help='Amount of RAM (GB) allocated to the Java heap. (2)')
    parser_b2f.add_argument('--PicardJar', '-jar', metavar='<path>',
        type=str, default=os.path.join(root_dir, 'libs/picard.jar'),
        help='Location of picard.jar (libs/picard.jar)')
    # Add function to call if selected
    parser_b2f.set_defaults(func=basecalls2fastq.run)
    
    # Create parser for the demultiplexer
    parser_demux = subparsers.add_parser('demux')
    # Required positional args
    parser_demux.add_argument('inDir', metavar='<inDir>',
        type=str, help='Directory created by bcl2fastq in output folder.')
    parser_demux.add_argument('sampleSheet', metavar='<SampleSheet.csv>',
        type=str, help='MiSeq SampleSheet.csv file, containing index info.')
    # Optional args
    parser_demux.add_argument('--uniqID', '-u', metavar='<str>', type=str,
        default=None, help='Unique ID to append to output folder. (None)')
    parser_demux.add_argument('--numCPU', '-p', metavar='<int>', type=int,
        default=1, help='Number of CPUs to use. (1)')
    parser_demux.add_argument('--minProb', '-min', metavar='<float>',
        type=float, default=0.05, help=('Minimum probability of a match else'
                                        ' discard. (0.05)'))
    parser_demux.add_argument('--phredOffset', '-s', metavar='<int>', type=int,
        required=False, default=33, help='FASTQ phred score offset (33)')
    parser_demux.add_argument('--indexQual', '-i', metavar='<int>', type=int,
        default=30, help='Phred-score given to barcode indexes (30)')
    # Add function to call if selected
    parser_demux.set_defaults(func=demultiplexer.run)

    # Add version number to the parser
    parser.add_argument('-v', '--version', action='version', version='v0.9.0')

    # Parse the arguments
    args = parser.parse_args()

    # Workaround for sub-parser bug (http://bugs.python.org/issue16308)
    try:
        a = getattr(args, "func")
    except AttributeError:
        print(parser.print_help())
        sys.exit(0)

    # Parse the arguments
    return args

if __name__ == '__main__':
    main()
