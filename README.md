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
$ bowtie -p 4 -S -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R1.fastq > input_1.sam
$ bowtie -p 4 -S -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 R2.fastq > input_2.sam

# sort sam files
$ samtools view -@ 4 -bS input_1.sam | samtools sort -n -@ 4 - sorted
$ samtools view -@ 4 sorted.bam > input_1.sam
$ samtools view -@ 4 -bS input_2.sam | samtools sort -n -@ 4 - sorted
$ samtools view -@ 4 sorted.bam > input_2.sam
$ rm sorted.bam

# analysis
$ cd /path/to/CASTIN
$ java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interactome.Main -p /path/to/CASTIN/inputdir/input -o /path/to/CASTIN/outputdir
```

### variable length single-end input (e.g., iontorrent)
```bash
# input preparation
$ cd /path/to/CASTIN/inputdir
$ bowtie -S --sam-nohead -a -v 1 /path/to/CASTIN/parameters/hg38_mm10/hg38_mm10 input.fastq > input.sam

# analysis
$ cd /path/to/CASTIN
$ java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interactome.Main -s /path/to/CASTIN/inputdir/input.sam -o /path/to/CASTIN/outputdir
```

## Options

```bash
-s single-end input file prefix (cannot be specified with -p)

-p paired-end input file prefix (cannot be specified with -s)

-o output directory

-d directionality for paired end input (0: undirectional 1: (forward, reversed) only, 2: (reversed, forward) only)

```

## Input files

CASTIN can take alignment sam files as input. Any alignment software may be used to produce the sam files, but bowtie aligner is recommended. Indexed reference sequence files for bowtie are included in the software package.

In case of using paired-end input, sort sam files by read names before processing.

## Output files

CASTIN produces multiple output files. All output files are written in the directory specified by the command line option “-o /path/to/CASTIN/outputdir”. 

### Interactive profile file 

- KEGGHPRD_result.txt :  This file contains ligand dependency (ligand ratio), receptor dependency (receptor ratio), and signal strength (average) for each interaction. It also contains the type of interaction and the link to the original pathway in KEGG. 

### Gene expression files

- Symbol_cancer.txt : Gene expression levels in cancer (human). The last column is the normalized gene expression level.
- Symbol_stroma.txt : Gene expression levels in stroma (mouse). The last column is the normalized gene expression level.
- Refseq_cancer.txt : Read count, the number of mappable reads and coverage information of each Refseq gene in cancer (human).
- Refseq_stroma.txt : Read count, the number of mappable reads and coverage information of each Refseq gene in storm (mouse).

### Image files

The following image files show the scatter plot before and after Poly-A or GC content derived bias.
Files contains _estim_ indicates that only genes used for bias parameter estimation were plot.

- from-poly-A_after_correction.png
- from-poly-A_before_correction.png
- from-poly-A_estim_after_correction.png
- from-poly-A_estim_before_correction.png
- GC_after_correction.png
- GC_before_correction.png
- GC_estim_after_correction.png
- GC_estim_before_correction.png

## License

CASTIN is released under the GNU General Public License (GPL).

