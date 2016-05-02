#! /bin/sh

# sample script for preparing sequence mapping using vmatch.
# execute this script in the *root* directory, after executing 'prepare_mapping.sh'

TAX_1=hg38
TAX_2=mm10

THREADS=2
MIN_SELF_MATCH_LENGTH=50

cd parameters/${TAX_1}_${TAX_2}

mkvtree -db ${TAX_1}_${TAX_2}_refMrna.fasta -dna -pl 12 -allout -v

vmatch -l 50 -v -d -p ${TAX_1}_${TAX_2}_refMrna.fasta > ${TAX_1}_${TAX_2}_selfmatch.txt

ruby ../../scripts/summarize_vmatch_result.rb ./${TAX_1}_${TAX_2}_refMrna.fasta ../refLink.txt ./${TAX_1}_${TAX_2}_selfmatch.txt dp > ./${TAX_1}_${TAX_2}_selfmatch_summary.txt
ruby ../../scripts/summarize_vmatch_result.rb ./${TAX_1}_${TAX_2}_refMrna.fasta ../refLink.txt ./${TAX_1}_${TAX_2}_selfmatch.txt d > ./${TAX_1}_${TAX_2}_selfmatch_summary_direct.txt

rm ./${TAX_1}_${TAX_2}_selfmatch.txt
