# MRC: Multiple reads set Clustering and Compression

Multiple Reads Sets Clustering (MRC) is a novel clustering approach for the compression of multiple FASTQ files which construct features space based on the freuqency of minimizers in each file.

## Download & Install

	git clone git@github.com:ttan6729/MRC.git
	cd MRC
	sh install.sh

## Usage
Usage:
Compression - compresses FASTQ datasets. Output written to '*.MRC' file
```
./MRC.sh -a m -r file.txt (compress with minicom, file contain name of to be compressed files)
./MRC.sh -a p -r file.txt (compress with PgRc, file contain name of to be compressed files)
```
Decompression - decompress compressed files, decompressed fastq files are written to '._MRC' folder
```
./MRC.sh -d file.MRC
```
Options:
        -r      compression mode
        -h              print help message
        -t              number of threads, default: 12
        -k              length of k-mer, k <= 10, default: 8
        -e              threshold percentage, default: 2.0
        -d              a compressed file .MRC [only for decompression]
        -t              number of threads, default: 12
        
##Example
```
./MRC.sh -a p -r test.txt
./MRc.sh -d testp.MRC
```
