#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author: Edward Mountjoy
#

import sys
import re
import os
from shutil import rmtree
import subprocess

def run(args):

    # Parse runParameters.xml file to get machine, flow, run names
    runParam = parse_runParameters(args.runParamXML)
    machineName, runBarcode, flowcellBarcode, readStructure = runParam

    # Check that top output folder exists
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    # Create unique output folder
    uniq_name = '{0}_{1}_{2}'.format(machineName, runBarcode, flowcellBarcode)
    out_dir = os.path.join(args.outDir, uniq_name)
    if os.path.exists(out_dir):
        response = input('{0} exists. Would you like to overwrite it? [y/n] '.format(out_dir))
        if response == 'y':
            rmtree(out_dir)
        else:
            sys.exit()
    os.makedirs(out_dir)
    # Create multiplexed folder inside unique output folder
    out_dir = os.path.join(out_dir, "multiplexed")
    os.makedirs(out_dir)

    # Make outprefix
    out_prefix = os.path.join(out_dir, 'multiplexed')

    # Build command
    cmd = ['java',
           '-Xmx{0}g'.format(args.JavaRAM),
           '-jar',
            args.PicardJar,
           'IlluminaBasecallsToFastq',
           'BASECALLS_DIR={0}'.format(args.baseCallDir),
           'LANE={0}'.format(args.lane),
           'OUTPUT_PREFIX={0}'.format(out_prefix),
           'RUN_BARCODE={0}'.format(runBarcode),
           'MACHINE_NAME={0}'.format(machineName),
           'FLOWCELL_BARCODE={0}'.format(flowcellBarcode),
           'READ_STRUCTURE={0}'.format(readStructure),
           'NUM_PROCESSORS={0}'.format(args.numCPU),
           'MAX_READS_IN_RAM_PER_TILE={0}'.format(args.readsPerTile),
           'MAX_RECORDS_IN_RAM={0}'.format(args.MaxInRam),
           #'COMPRESS_OUTPUTS=true',
           
           # The below options should be uncommented for debugging as we do
           # not need all reads to be extracted. Compressed files do not work
           # if a tile limit is used
           'COMPRESS_OUTPUTS=false',
           'TILE_LIMIT=1' ,
           #'VERBOSITY=DEBUG'
          ]

    # Run IlluminaBasecallsToFastq
    subprocess.call(cmd)

    return 0

def parse_runParameters(filen):
    """ Need to get: MachineName = ScannerID
                     RunBarcode = RunNumber
                     FlowcellBarcode = Barcode
                     ReadStructure
    """
    # Compile regex patterns
    re_patterns = {'MachineName':re.compile(r'<ScannerID>(.+)</ScannerID>'),
                   'RunBarcode':re.compile(r'<RunNumber>(.+)</RunNumber>'),
                   'FlowcellBarcode':re.compile(r'<Barcode>(.+)</Barcode>')}
    re_struct = re.compile(r'<RunInfoRead Number="([0-9]+)" NumCycles="([0-9]+)" IsIndexedRead="([Y|N])" />')
    labels = {}
    structure = {}

    # Iterate over each line and check for matches
    with open(filen, 'r') as in_handle:
        for line in in_handle:
            # Check for machine, run and cell barcodes
            for label in re_patterns:
                match_obj = re_patterns[label].search(line)
                # If a match is found
                if match_obj:
                    labels[label] = match_obj.group(1)
            # Check for read structure
            match_obj = re_struct.search(line)
            if match_obj:
                num = int(match_obj.group(1))
                cycles = int(match_obj.group(2))
                index = match_obj.group(3)
                structure[num] = (cycles, index)

    # Check that each was found
    for label in re_patterns:
        if label not in labels:
            sys.exit('Could not find {0} in {1}'.format(label, filen))

    # Build read structure
    rs = []
    for key in sorted(structure.keys()):
        if structure[key][1] == 'Y':
            rtype = "B"
        elif structure[key][1] == 'N':
            rtype = "T"
        rs.append("{0}{1}".format(structure[key][0], rtype))
    rs = "".join(rs)

    return (labels['MachineName'], labels['RunBarcode'],
            labels['FlowcellBarcode'], rs)
