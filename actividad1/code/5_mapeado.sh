#!/usr/bin/env bash

## Mapeado
if [ ! -f ../results/mapeado/SRR13795616_hisat2.sam ]
then
   ## Primero indexamos el genoma de referencia de la
   ## especie de estudio: Mus musculus
   hisat2-build  ../data/raw/genoma.fa ../data/reference_genome/genome

   ## Realizamos el alineamiento con el genoma de referencia indexado
   hisat2 -k2 -U ../data/processed/ejercicio2_trimmed/SRR13795616_trimmed.fq.gz -x ../data/reference_genome/genome --rna-strandness R -S ../results/mapeado/SRR13795616_hisat2.sam
else
   echo "Los datos ya estaban mapeados"
fi

