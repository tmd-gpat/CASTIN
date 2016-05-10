# CASTIN: CAncer-STromal INteractome analysis

[![DOI](https://zenodo.org/badge/21913/tmd-gpat/CASTIN.svg)](https://zenodo.org/badge/latestdoi/21913/tmd-gpat/CASTIN)

## General Notes

CASTIN (CAncer-STromal INteractome analysis) is a framework for the evaluation of cancer-stromal interactome from RNA-seq data using cancer xenograft models.
For each ligand-receptor interaction which is derived from curated protein-protein interaction database, CASTIN summarizes gene expression profiles of tumor and stroma into three evaluation indices.

## Requirements

CASTIN is depending on the following software/libraries.

- JDK (>= 1.8)
- Ruby (>= 2.0)
- R (>= 3.1.2)
- rJava

Following software are needed for preparing input and parameter files.

- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [vmatch](http://www.vmatch.de/)

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
$ export JRI_DIR=/path/to/R/library/rJava/jri/
$ export R_HOME=/path/to/R/home/
```

## Basic Usage

### 100bp paired-end input (e.g., illumina HiSeq)
```bash
# input preparation
$ cd /path/to/CASTIN/inputdir
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R1.fastq > input_100_1.sam
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R2.fastq > input_100_2.sam

# analysis
$ cd /path/to/CASTIN
$ java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interaction.Main -p /path/to/CASTIN/inputdir/input -l 100 -o /path/to/CASTIN/outputdir
```

### variable length single-end input (e.g., iontorrent)
```bash
# input preparation
$ cd /path/to/CASTIN/inputdir
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 input.fastq > input.sam

# analysis
$ cd /path/to/CASTIN
$ java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interaction.Main -s /path/to/CASTIN/inputdir/input.sam -o /path/to/CASTIN/outputdir
```

## License

CASTIN is released under the GNU General Public License (GPL).

