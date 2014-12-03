aware demultiplexer
===================

Current version: v0.9.1

Probabilistic demultiplexer for Illumina bcl files. Works with single or dual-
indexed reads, and single or pair-end reads.

## Setup

#### Dependancies
- Python (>=3.2)
- Java Runtime Environment (tested with openjdk-7-jre and Java SE Runtime Env 1.8)

#### Download
The repository can be downloaded using git `git clone https://github.com/edm1/aware-demultiplexer.git` or by following the *Download ZIP* link on the right.

#### Recommended
- [PyPy3 2.4.0](http://pypy.org/)

Using pypy3 instead of python3 will give approximately 3x speed up.

## Usage

The script is split into two sub-commands `bcl2fastq` and `demux`.

```
usage: pypy3 aware.py [-h] [-v] <subcommand> [options]

The aware.py sub-commands include:
    bcl2fastq    Extracts multiplexed reads and barcodes from Illumina bcl
                 files.
    demux        Demultiplex the fastqs extracted by bcl2fastq using indexes
                 provided in sampleSheet.csv.
```

To see further help for each sub-command use `pypy3 aware.py <subcommand> -h`.

#### Sub-command: bcl2fastq

`bcl2fastq` extracts read and barcode fastq files from Illumina .bcl files.

```
usage: pypy3 aware.py bcl2fastq [-h] [--outDir <str>] [--numCPU <int>]
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
If you get the error "Could not find a format with available files for the following data types: Position", it is because the folder `MiSeqOutput/InterOp` needs to be present in addition to `MiSeqOutput/Data/Intensities/BaseCalls`.

## Step 2: Demultiplex FASTQs into separate sample files
Barcodes are checked against sample indexes. The probability that the two underlaying sequences (barcode and index) are the same is calculated using the phred quality scores. I.e. the probability that the true bases match given that the sequenced bases match/mismatch is calculated across the index/barcode. Script is run using:

```
python [or pypy] 2_demultiplex_fastqs.py [-h] --InputTempDir <dir> --SampleSheet
                                              <SampleSheet.csv> --OutID <str>
                                              [--ResultsDir <dir>] [--MinProb <float>]
                                              [--PhredOffset <int>] [--IndexQual <int>]

```

##### Required
- `--InputTempDir` - Directory containing temp files output by 1_run_IlluminaBasecallsToFastq.py
- `--SampleSheet` - MiSeq output SampleSheet.csv file containing barcode indexes and sample names.
- `--OutID` - Unique ID to append to the output folder name.

##### Optional
- `--MinProb` - Minimum probability of a match, else discard. (0.05)
- `--IndexQual` - Phred-score given to index sequence. (30)
- `--OutDir` - Results directory. (default: ./results)
- `--PhredOffset` - FASTQ phred score offset. (33)

### Output
For each of the sample indexes + not-assigned reads, the script outputs:
1. FASTQs for the demultiplexed reads
2. FASTQs containing the barcodes of each of the reads
