#! /bin/sh

# execute this script on the *root* directory!
# automatic script for downloading hg38(cancer) and mm10(stroma) sequence data
#
# [updated at Apr.29 2016]

# make parameter directories
if [ ! -d parameters ]; then
    mkdir parameters
fi
if [ ! -d parameters/hg38 ]; then
    mkdir parameters/hg38
fi
if [ ! -d parameters/mm10 ]; then
    mkdir parameters/mm10
fi

# refLink
cd parameters
wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refLink.txt.gz"
gunzip refLink.txt.gz

# homologene.data
wget "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"

# refMrna / refNames / refSeqLen
cd hg38
wget "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz"
gunzip refMrna.fa.gz
ruby ../../scripts/refMrna2refNames.rb < refMrna.fa > refNames.txt
ruby ../../scripts/refMrna2refSeqLen.rb < refMrna.fa > refSeqLen.txt
cd ..

cd mm10
wget "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/refMrna.fa.gz"
gunzip refMrna.fa.gz
ruby ../../scripts/refMrna2refNames.rb < refMrna.fa > refNames.txt
ruby ../../scripts/refMrna2refSeqLen.rb < refMrna.fa > refSeqLen.txt
cd ..
