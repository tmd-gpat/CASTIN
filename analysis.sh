#! /bin/sh
# use this sample script for by replacing paths.

# supposed directory structure (paired-end)
#   CASTIN/ (cwd)
#     bin/
#     lib/
#     parameters/
#     ...
#
#   /path/to/CASTIN/input/dir/
#     input_prefix_1_100.sam
#     input_prefix_2_100.sam
#   /path/to/CASTIN/output/dir/
#
# * input sam files should be named as "PREFIX_1_LENGTH.sam" and "PREFIX_2_LENGTH.sam".

# set the path of your JRI (rJava) install directory
JRI_DIR=/usr/lib64/R/library/rJava/jri/

java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interaction.Main -p /path/to/CASTIN/input/dir/input_prefix -l 100 -o /path/to/CASTIN/output/dir/

