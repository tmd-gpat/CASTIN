# CASTIN: CAncer-STromal INteractome analysis

## General Notes

## Requirements

CASTIN is depending on the following software/libraries.

- JDK (>= 1.8)
- Ruby (>= 2.0)
- R (>= 3.1.2)
- rJava

Following software are needed for preparing input and parameter files.

- bowtie
- vmatch

## Installation

```bash
$ git clone git@github.com:tmd-gpat/CASTIN.git
$ cd CASTIN
$ ln -s /path/to/your/JRI.jar lib/JRI.jar
$ ant
```

Set appropriate *JAVA_HOME* variable beforehand.

## Preparation

```bash
$ cd CASTIN
$ ./scripts/download_parameter_files.sh # download reference sequences, etc.
$ ./scripts/prepare_mapping.sh          # make combined reference sequence and bowtie index
$ ./scripts/calculate_selfmatching.sh   # prepare self-matching list for mappability correction
```

## Basic Usage

### 100bp paired-end input
```bash
# input preparation
$ cd /path/to/input/dir
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R1.fastq > input_100_1.sam
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R2.fastq > input_100_2.sam

# analysis
$ cd /path/to/CASTIN
$ export JRI_DIR=/path/to/R/library/rJava/jri/
$ java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interaction.Main -p /path/to/input/dir/input -l 100 -o /path/to/output/dir
```

## Acknowledgement

