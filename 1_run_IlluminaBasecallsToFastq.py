#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# python run_IlluminaBasecallsToFastq.py --BaseCallDir ../test_data/MiSeqOutput/140321_M01520_0085_000000000-A7FD0/Data/Intensities/BaseCalls/ --RunParamXML ../test_data/MiSeqOutput/140321_M01520_0085_000000000-A7FD0/runParameters.xml --Lane 1 --ReadStructure 151T8B8B
#

import argparse
import re
import sys
import os
import inspect
from shutil import rmtree
import subprocess

def main():

    # Get root of project directory
    global root_dir
    root_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    #root_dir = os.path.split(cur_dir)[0]

    # Parse args
    args = parse_arguments()

    # Parse runParameters.xml file to get machine, flow, run names
    machineName, runBarcode, flowcellBarcode = parse_runParameters(args.RunParamXML)

    # Make an ouput directory
    out_dir = os.path.join(root_dir, 'temp')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_dir = os.path.join(out_dir, '{0}_{1}_{2}'.format(machineName,
                                                             runBarcode,
                                                             flowcellBarcode))
    if os.path.exists(out_dir):
        response = raw_input('{0} exists. Would you like to overwrite it? [y/n] '.format(out_dir))
        if response == 'y':
            rmtree(out_dir)
        else:
            sys.exit()
    os.makedirs(out_dir)

    # Make outprefix
    out_prefix = os.path.join(out_dir, 'raw_extract')

    # Build command
    cmd = ['java -Xmx{0}g -jar'.format(args.JavaRAM),
            args.JarLoc,
           'IlluminaBasecallsToFastq',
           'BASECALLS_DIR={0}'.format(args.BaseCallDir),
           'LANE={0}'.format(args.Lane),
           'OUTPUT_PREFIX={0}'.format(out_prefix),
           'RUN_BARCODE={0}'.format(runBarcode),
           'MACHINE_NAME={0}'.format(machineName),
           'FLOWCELL_BARCODE={0}'.format(flowcellBarcode),
           'READ_STRUCTURE={0}'.format(args.ReadStructure),
           'NUM_PROCESSORS={0}'.format(args.numCPU),
           'MAX_READS_IN_RAM_PER_TILE={0}'.format(args.readsPerTile),
           'MAX_RECORDS_IN_RAM={0}'.format(args.MaxInRam),
           'COMPRESS_OUTPUTS=true',
           
           # The below options should be uncommented for debugging as we do
           # not need all reads to be extracted. Compressed files do not work
           # if a tile limit is used
           #'COMPRESS_OUTPUTS=false',
           #'TILE_LIMIT=1' ,
           #'VERBOSITY=DEBUG'
          ]
    cmd = ' '.join(cmd)

    # Run IlluminaBasecallsToFastq
    subprocess.call(cmd, shell=True)

def parse_runParameters(filen):
    """ Need to get: MachineName = ScannerID
                     RunBarcode = RunNumber
                     FlowcellBarcode = Barcode
    """
    # Compile regex patterns
    re_patterns = {'MachineName':re.compile(r'<ScannerID>(.+)</ScannerID>'),
                   'RunBarcode':re.compile(r'<RunNumber>(.+)</RunNumber>'),
                   'FlowcellBarcode':re.compile(r'<Barcode>(.+)</Barcode>')}
    matches = {}

    # Iterate over each line and check for matches
    with open(filen, 'r') as in_handle:
        for line in in_handle:
            for label in re_patterns:
                match_obj = re_patterns[label].search(line)
                # If a match is found
                if match_obj:
                    matches[label] = match_obj.group(1)

    # Check that each was found
    for label in re_patterns:
        if label not in matches:
            sys.exit('Could not find {0} in {1}'.format(label, filen))

    return matches['MachineName'], matches['RunBarcode'], matches['FlowcellBarcode']


def parse_arguments():
    """ Load the arguments required to run IlluminaBasecallsToFastq.
    """

    parser = argparse.ArgumentParser(description='Run IlluminaBasecallsToFastq using picard.')

    # Required args
    parser.add_argument('--BaseCallDir',
                        metavar='<dir>',
                        type=str,
                        required=True,
                        help='Directory containing base call intensitites')
    parser.add_argument('--RunParamXML',
                        metavar='<runParameters.xml>',
                        type=str,
                        required=True,
                        help='runParameters.xml file')
    parser.add_argument('--Lane',
                        metavar='<int>',
                        type=int,
                        required=True,
                        help='Lane number')
    parser.add_argument('--ReadStructure',
                        metavar='<str>',
                        type=str,
                        required=True,
                        help='Description of the logical structure of clusters in an Illumina Run, e.g. 151T8B8B')

    # Optional arguments
    parser.add_argument('--JarLoc',
                        metavar='<*.jar>',
                        type=str,
                        required=False,
                        default=os.path.join(root_dir, 'lib/picard-tools-1.125/picard.jar'),
                        help='Location of picard.jar (../lib/picard-tools-1.125/picard.jar)')
    parser.add_argument('--numCPU',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=1,
                        help='Number of CPUs to use. (1)')
    parser.add_argument('--readsPerTile',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=120000,
                        help='Max number of reads in RAM per tile, reduce if you have problems with memory. (120000)')
    parser.add_argument('--MaxInRam',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=500000,
                        help='Maximum number of records that are stored in the RAM. (500000)')
    parser.add_argument('--JavaRAM',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help='Amount of RAM allocated to the Java heap. (2)')
    parser.add_argument('--version', action='version', version='v0.2')


    return parser.parse_args()

if __name__ == '__main__':
    main()
