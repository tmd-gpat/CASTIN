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
$ ./scripts/download_parameter_files.sh
$ ./scripts/prepare_mapping.sh
$ ./scripts/calculate_selfmatching.sh
```

## Basic Usage

## Acknowledgement

