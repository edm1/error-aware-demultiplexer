#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Edward Mountjoy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import sys
import argparse
import os
from sys import platform as _platform
import subprocess
from shutil import copyfile

def main():
    """ Installs pypy3.
    """

    # Parse the command line args
    args = parse_arguments()

    print("Installing...")
    install_pypy(args)

    return 0

def install_pypy(args):
    """ Function get and install pypy3 binary.
    """
    
    # Make input python 2.7 compatible
    if sys.version_info[0] >= 3:
        get_input = input
    else:
        get_input = raw_input

    # Confirm paths
    exit_msg = "\nExiting. Use --help to view install options."
    for msg in ["> Install path: {0} [y/n] ".format(args.dir),
                "> Bashrc path: {0} [y/n] ".format(args.bashrc)]:
        ret = get_input(msg)
        if not ret == "y":
            sys.exit(exit_msg)

    # Make output folder
    make_folders(args.dir)

    # Get and extract pypy3
    temp_pypy = "pypy3_2.4.0_download.tar.bz2"
    cmd = []
    if _platform == "linux" or _platform == "linux2":
        url = "https://bitbucket.org/pypy/pypy/downloads/pypy3-2.4.0-linux64.tar.bz2"
        cmd.append('wget {0} -O {1}'.format(url, temp_pypy))
    elif _platform == "darwin":
        url = "https://bitbucket.org/pypy/pypy/downloads/pypy3-2.4.0-osx64.tar.bz2"
        # OS X
        cmd.append('curl -o {0} -L {1}'.format(temp_pypy, url))
    # Unzip file
    cmd.append('tar -jxvf {0} --strip 1 -C {1}'.format(temp_pypy, args.dir))
    # Run command
    ret = subprocess.call(";".join(cmd), shell=True)
    if not ret == 0:
        sys.exit("There was a problem downloading or extracting pypy. Exiting.")
    # Remove download
    os.remove(temp_pypy)

    # Create backup of bashrc
    bashrc_backup = "{0}_backup".format(args.bashrc)
    if os.path.exists(args.bashrc):
        copyfile(args.bashrc, bashrc_backup)
        print("\nCreated backup for of {0} at {1}.".format(args.bashrc, bashrc_backup))
    # Add pypy3 bin to PATH
    pypy_bin = os.path.join(args.dir, "bin")
    lines = ["\n# PyPy3 2.4.0 bin PATH - created by aware-demultiplexer",
             "export PATH=$PATH:{0}\n".format(pypy_bin)]
    with open(args.bashrc, 'a') as out_h:
        for line in lines:
            out_h.write(line + "\n")

    print("Finished installing PyPy3")

def make_folders(outDir):
    # Get list of folders that need checking
    check_dirs = []
    check_dir = outDir
    while True: #not :
        # Check that its not home dir
        try:
            if os.path.samefile(check_dir, os.getenv("HOME")):
                break
        except FileNotFoundError:
            pass
        # Append file
        check_dirs.append(check_dir)
        check_dir = os.path.split(check_dir)[0]
    # Check those folders
    for check_dir in check_dirs[::-1]:
        if not os.path.exists(check_dir):
            os.makedirs(check_dir)
    return 0

def parse_arguments():
    """ Will parse the command line arguments arnd return the arg object.
    """
    home_dir = os.getenv("HOME")
    parser = argparse.ArgumentParser(
        description="Installs PyPy3 2.4.0 in user's home directory")
    parser.add_argument("--dir", metavar='<installDir>',
        help="Directory to install PyPy3 to. (Default: ~/programs/pypy3-2.4.0)",
        default=os.path.join(*[home_dir, "programs", "pypy3-2.4.0"]))
    parser.add_argument("--bashrc", metavar='<bashrc>',
        help=("Location of basrc file (or equivalent) to append pypy3 bin path "
              "to. (Default: ~/.bashrc)"),
        default=os.path.join(home_dir, ".bashrc"))

    # Parse the arguments
    return parser.parse_args()

if __name__ == '__main__':
    main()
