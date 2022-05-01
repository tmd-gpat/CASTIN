#!/bin/bash

echo "create output directory $2"
mkdir /data/$2

java -cp "./bin:./lib/*" -Xmx16g -Xms8g -Djava.library.path=$JRI_DIR interactome.Main -p /data/$1 -o /data/$2
