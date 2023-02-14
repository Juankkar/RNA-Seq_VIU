#!/usr/bin/env bash

## Descargamos el archivo de las notaciones en caso de no ternerlo
if [ ! -f ../data/annotation/genome.gtf ]
then 
    wget -P ../data/annotation https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
    gunzip ../data/annotation/Mus_musculus.GRCm38.102.gtf.gz
    mv ../data/annotation/Mus_musculus.GRCm38.102.gtf ../data/annotation/genome.gtf
    echo "Se han descargado los datos correctamente"
else
    echo "Los datos ya estaban descargados"
fi

## Y en caso de no tener la matriz de recuentos ya creada lo hacemos
if [ ! -f ../results/tables/SRR13795616_counts.tsv ]
then
    htseq-count -f bam -r pos -m union -s reverse -t exon -i gene_id ../results/mapeado/SRR13795616_hisat2.sorted.bam ../annotation/genome.gtf > ../results/tables/SRR13795616_counts.tsv
    echo "Se ha creado una matriz de recuentos"
else
    echo "La matriz de recuentos ya estaba creada"
fi
