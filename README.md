Error Aware Demultiplexer (EAD)
===================

Current version: v1.0.3

EAD is a probabilistic demultiplexer for Illumina BCL files. It works with:
- single or dual-indexed reads
- single or pair-end reads

The demultiplexer has only been tested with MiSeq BCLs but theoretically should work for any Illumina platform.

## Setup

#### Dependancies
- Python (>=3.2)
- Java Runtime Environment (tested with openjdk-7-jre and Java SE Runtime Env 1.8)

#### Download
The repository can be downloaded using `git clone git@github.com:edm1/error-aware-demultiplexer.git` or by following the *Download ZIP* link on the right.

#### Recommended
- [PyPy3](http://pypy.org/) (>=2.4.0)

Using pypy3 instead of python3 will give approximately 3x speed up. If you do not have pypy3 installed, it can be installed using `python extras/install_pypy3-2.4.0.py`.

## Usage

The script is split into two sub-commands `bcl2fastq` and `demux`.

```
usage: pypy3 ead.py [-h] [-v] <subcommand> [options]

The ead.py sub-commands include:
    bcl2fastq    Extracts multiplexed reads and barcodes from Illumina bcl
                 files.
    demux        Demultiplex the fastqs extracted by bcl2fastq using indexes
                 provided in sampleSheet.csv.
```

To see further help for each sub-command use `pypy3 ead.py <subcommand> -h`.

#### Sub-command: bcl2fastq

`bcl2fastq` extracts read and barcode fastq files from Illumina .bcl files.

```
usage: pypy3 ead.py bcl2fastq [-h] [--outDir <str>] [--numCPU <int>]
                                [--readsPerTile <int>] [--MaxInRam <int>]
                                [--JavaRAM <int>] [--PicardJar <path>]
                                <baseCallDir> <runParameters.xml> <lane>
```


##### Required
- `<baseCallDir>` - `MiSeqOutput/Data/Intensities/BaseCalls` directory containing BCL files.
- `<runParameters.xml>` - `MiSeqOutput/runParameters.xml` file.
- `<lane>` - Lane number

##### Optional
- `--outDir` - Location to create output files. (./output)
- `--numCPU` - Number of CPUs to use. (1)
- `--readsPerTile` - Max number of reads in RAM per tile, reduce if you have problems with memory. (120000)
- `--MaxInRam` - Maximum number of records that are stored in the RAM. (500000)
- `--JavaRAM` - Amount of RAM allocated to Java heap. Increase if having problems. (2)
- `--JarLoc` - Location of picard.jar (./libs/picard.jar)

##### Troubleshooting
If you get the error "Could not find a format with available files for the following data types: Position", it is because the folder `MiSeqAnalysis/InterOp` needs to be present in addition to `MiSeqAnalysis/Data/Intensities/BaseCalls`. I.e. the following starred (*) folders/files are required:

```
MiSeqAnalysis
├── Config
├── Data
│   ├── Intensities
│   │   ├──  BaseCalls -- *
│   │   ├──  L001 ------- *
│   │   ├── [L002] ------ *
│   │   ├── [L...] ------ *
│   │   └──  Offsets
│   ├── RTA Logs
│   └── Tile Status
├── Images
├── InterOp ------------ *
├── Logs
├── Recipe
├── Thumbnail_Images
├── RunInfo.xml
├── runParameters.xml -- *
└── SampleSheet.csv ---- *
```

#### Sub-command: demux

Demultiplex the fastqs extracted by bcl2fastq using index sequences provided in sampleSheet.csv.

```
usage: pypy3 ead.py demux [-h] [--uniqID <str>] [--minProb <float>]
                            [--phredOffset <int>] [--indexQual <int>]
                            <inDir> <SampleSheet.csv>
```

##### Required
- `<inDir>` - Directory containing `multiplexed` folder output by bcl2fastq.
- `<sampleSheet.csv>` - MiSeq SampleSheet.csv file containing barcode indexes and sample names.

##### Optional
- `--uniqID` - Unique ID to append to output folder. Useful if testing multiple parameters. (None)
- `--minProb` - Minimum overall probability that a barcode and index match, else sample = "not_assigned". (0.05)
- `--phredOffset` - FASTQ phred score offset. (33)
- `--indexQual` - Phred-score given to all bases in index sequence. (30)

## Acknowledgements

The Broad Institutes [Picard tool](github.com/broadinstitute/picard) is used to extract FASTQs from Illumina BCL files. I have redistributed it unchanged in `./libs/picard.jar`.

## License

```
The MIT License (MIT)

Copyright (c) 2014 Edward Mountjoy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

The function src.fastqparser.fastqIterator was adpated from the biopython source.

```
Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
```
