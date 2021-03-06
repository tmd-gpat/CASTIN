#! /bin/sh

# sample script for preparing sequence mapping using bowtie.
# execute this script in the *root* directory.

TAX_1=hg38
TAX_2=mm10

# concat reference sequences
cd parameters
mkdir ${TAX_1}_${TAX_2}
cat ${TAX_1}/refMrna.fa ${TAX_2}/refMrna.fa > ${TAX_1}_${TAX_2}/${TAX_1}_${TAX_2}_refMrna.fasta

# build bowtie and bowtie2 indexes
cd ${TAX_1}_${TAX_2}
bowtie-build -f ${TAX_1}_${TAX_2}_refMrna.fasta ${TAX_1}_${TAX_2}

