#! /bin/sh

# sample script for preparing sequence mapping using bowtie.
# execute this script in the *root* directory, after 'prepare_mapping.sh'

TAX_1=hg38
TAX_2=mm10

THREADS=2
MIN_SELF_MATCH_LENGTH=50

cd parameters/${TAX_1}_${TAX_2}
bowtie2-build ${TAX_1}_${TAX_2}_refMrna.fasta ${TAX_1}_${TAX_2}

bowtie2 -p $THREADS -x ${TAX_1}_${TAX_2} --local -a --score-min L,${MIN_SELF_MATCH_LENGTH},0 --ma 1 --mp 10000,10000 --rdg 10000,10000 --rfg 10000,10000 --no-hd -f -U ${TAX_1}_${TAX_2}_refMrna.fasta | awk '{if ($1!=$3) print $1, $3, $4, $6, $19 }' > ${TAX_1}_${TAX_2}_selfmatching.txt

