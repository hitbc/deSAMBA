deSAMBA
======

deSAMBA: fast and accurate classification of metagenomics long reads with sparse approximate matches

## Table of Contents
1. [Dependency](#dependency)
2. [Quick start](#Quick-start)
3. [Introduction](#Introduction)
4. [Memory usage](#Memory-usage)
5. [Build project](#build-project)
6. [Build index](#build-index)
7. [Run classifation](#run-classifation)
8. [Run analysis](#run-analysis)
9. [Demo datat](#Demo-data)

## Dependency

Jelly fish depends on g++, make, automake and libtool; deSAMBA depend on zlib.

Run following commands to install them (Ubuntu).
```
sudo apt-get install zlib1g-dev
sudo apt-get install automake
sudo apt-get install libtool
sudo apt-get install make
sudo apt-get install g++
```
## Quick start
```
git clone https://github.com/hitbc/deSAMBA.git --depth=1
cd ./deSAMBA
bash ./build
bash ./build-index ./demo/viral-gs.fa ./demo_index
#classify
./bin/deSAMBA classify -t 4 ./demo_index ./demo/ERR1050068.fastq -o ./ERR1050068.sam
#analysis
./bin/deSAMBA analysis ana_meta ./ERR1050068.sam ./demo/nodes.dm
```

## Introduction

DeSAMBA(de Bruijn graph-based Sparse Approximate Match Block Analyzer) is a tailored long read classification approach using a novel sparse approximate match-based pseudo alignment algorithm. Benchmarks on real datasets demonstrate that deSAMBA enables to simultaneously achieve fast speed and good classification yields, which outperforms state-of-the-art tools and has many potentials to cutting-edge metagenomics studies. 

## Memory usage

DeSAMBA fits modern servers and the peak memory footprint depends on the size and complexity of reference genome. Memory used for building index is bigger than running classification, so you can build index in a server with big memory and distribute the index in a small one. For option "all"(see build index options part for detail), 160 Gigabytes memory is required to build index, and 69 Gigabytes to run classification (Nov 2018, 35 Gigabytes of reference genome). For option "viral", 3 Gigabytes memory is required to build index, and 1 Gigabytes to run classification. As NCBI RefSeq database growing fast, more memory mey be required at present. 

## Build project

To build Jellyfish(v1.10) and deSAMBA
```
git clone https://github.com/hitbc/deSAMBA.git --depth=1
cd PATH/deSAMBA
bash ./build
```
## Build index

**Usage**
  
    build-index <REF_DIR> <index_dir>
   
**Options**
```
   <ref_file>  FILE    "all""viral" OR [file] of reference sequences.
   <IDX_DIR>   FOLDER  folder to store index.
```
   all: Using lastest NCBI RefSeq bacteria+viral+archaea database.
                                
   viral: Using lastest NCBI RefSeq viral database.
                             
   [file]: Otherwise using user defined reference sequences[.fa], storing them in one file.
   When you have more than one files,	combined them into one.
   using command:
```
        find fasta_file_dir/ -name "*.fasta"   | xargs -n 1 cat > WGS_FILE.fa
```
   to combined mulity fasta files into one.
 
**Demo**

Run "build-index" to build index, one of following three commands can be chosen
```
build-index all ./index_dir
build-index viral ./index_dir
build-index [file].fa ./index_dir
```

## Run classifation
```
/bin/deSAMBA classify -t 4 ./index_dir read.fastq -o result.sam
```

```
  Usage:     deSAMBA  classify  [Options] <IndexDir> [ReadFiles.fa][...]
  Basic:   
    <IndexDir>      FOLDER   the directory contains deSAMBA index
    [ReadFiles.fa]  FILES    reads files, FASTQ(A) format, separated by space
  Options:
    -h,             help
    -t, INT         number of threads[4]
    -l, INT         minimum matching length(will be ignored when classifying short NGS reads)[170]
    -r, INT         max Output number of secondary alignments[5]
    -o, FILE        output results into file [stdout]
    -s, FILE        MIN score(the number of 9-mer match)[64]
    -f, STR         output format, one of:
                    - SAM: SAM without SEQ and QUAL, default
                    - SAM_FULL: normal SAM
                    - DES: smallest format
                    - DES_FULL: all results will be showed, ignore '-r' option

```
Increasing "-l" and "-s" options will increase accuracy but decrease sensitivity.
[read.fastq] can be long noisy reads(error rate < 25%) or short NGS reads.

## Run analysis
```
./download taxnomy #download node.dmp
./bin/deSAMBA analysis ana_meta ./result.sam ./nodes.dmp 
```
## Demo data

We use part of viral whole genome sequences as demo reference(NCBI viral reference sequences, AUG 2018). And ERR1050068.fastq 
(Zaire ebolavirus, https://www.ncbi.nlm.nih.gov/sra/ERR1050068, tid = 186538) as demo meta-genomic data set.
Besides, we provide nodes.dump (AUG 2018) to running analysis.
```
cd demo
cd ..
#build index
bash ./build-index ./demo/viral-gs.fa ./demo_index
#classify
./bin/deSAMBA classify -t 4 ./demo_index ./demo/ERR1050068.fastq -o ./ERR1050068.sam
#analysis
./bin/deSAMBA analysis ana_meta ./ERR1050068.sam ./demo/nodes.dmp
```


