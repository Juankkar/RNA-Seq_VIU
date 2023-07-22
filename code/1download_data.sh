#!/usr/bin/env bash

sample_list=$(cut -f $1 metadata/report.tsv | tail -n 8)

for sample in ${sample_list[*]}
do
    wget -P data/raw/ $sample
done
