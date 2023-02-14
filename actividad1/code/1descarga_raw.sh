#!/usr/bin/env bash

echo ""
echo "###############################"
echo "Descarga de los archivos crudos"
echo "###############################"
echo ""
## En caso de no tener los datos crudos los descargamos
if [ ! -f ../data/raw/SRR13795616.fastq.gz ] && [ ! -f ../data/raw/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz ]
then
    wget -P ../data/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/016/SRR13795616/SRR13795616.fastq.gz 
    wget -P ../data/raw ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
    gzip -d ../data/raw/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
    mv ../data/raw/Mus_musculus.GRCm38.dna_sm.toplevel.fa ../data/raw/genoma.fa
    echo "Los datos han sido descargados en el directorio ../data/raw/"
else
    echo "Los datos ya estaban descargados en ../data/raw/"
fi

