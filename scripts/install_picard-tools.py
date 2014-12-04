#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# python run_IlluminaBasecallsToFastq.py --BaseCallDir ../test_data/MiSeqOutput/140321_M01520_0085_000000000-A7FD0/Data/Intensities/BaseCalls/ --RunParamXML ../test_data/MiSeqOutput/140321_M01520_0085_000000000-A7FD0/runParameters.xml --Lane 1 --ReadStructure 151T8B8B
#

import os
import inspect
import subprocess
from sys import platform as _platform

def main():

    # Get root of project directory
    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    root_dir = os.path.split(cur_dir)[0]

    # Make lib folder if it doesn't exist
    out_dir = os.path.join(root_dir, 'lib')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    os.chdir(out_dir)

    # Build the command
    cmd = []
    url = 'https://github.com/broadinstitute/picard/releases/download/1.125/picard-tools-1.125.zip'
    filen = url.rsplit("/", 1)[1]
    if _platform == "linux" or _platform == "linux2":
        # linux
        cmd.append('wget {0}'.format(url))
    elif _platform == "darwin":
        # OS X
        cmd.append('curl -OL {0}'.format(url))
    # Unzip file
    cmd.append('unzip -o {0}'.format(filen))
    # Run command
    subprocess.call(";".join(cmd), shell=True)

    # Tidy up
    os.remove(filen)

if __name__ == '__main__':
    main()
