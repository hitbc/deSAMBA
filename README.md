deSAMBA-meta
======

## Table of Contents
1. [dependency](#dependency)
2. [build project](#build-project)
3. [build index](#build-index)
4. [run classifation](#run-classifation)
5. [run analysis](#run-analysis)
6. [Demo datat](#Demo-data)

## dependency

Jelly fish depends on g++, make, automake and libtool; deSAMBA depend on zlib.

Run following commands to install them (Ubuntu).
```
sudo apt-get install zlib1g-dev
sudo apt-get install automake
sudo apt-get install libtool
sudo apt-get install make
sudo apt-get install g++
```
## build project

To build Jellyfish(v1.10) and deSAMBA
```
git clone https://github.com/hitbc/deSAMBA.git --depth=1
cd PATH/deSAMBA
bash ./build
```
## build index

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
 
**demo**
Run "build-index" to build index, one of following three commands can be chosen
```
build-index all ./index_dir
build-index viral ./index_dir
build-index [file].fa ./index_dir
```

## run classifation
```
/bin/deSAMBA classify -t 4 ./index_dir read.fasta -o result.sam
```
## run analysis
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


