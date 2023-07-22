#!/usr/bin/env bash

wget -P data/annotations/ \
    -O data/annotations/genome.gtf.gz \
    https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz

gzip -d data/annotations/genome.gtf.gz

