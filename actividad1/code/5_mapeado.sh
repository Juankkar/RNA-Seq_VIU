#!/usr/bin/env bash

## Mapeado
if [ ! -f ../results/mapeado/SRR13795616_hisat2.sam ]
then
   ## Primero indexamos el genoma de referencia de la
   ## especie de estudio: Mus musculus
   hisat2-build  ../data/raw/genoma.fa ../data/reference_genome/genome

   ## Realizamos el alineamiento con el genoma de referencia indexado
   hisat2 -k1 -U ../data/processed/ejercicio2_trimmed/SRR13795616_trimmed.fq.gz -x ../data/reference_genome/genome --rna-strandness R -S ../results/mapeado/SRR13795616_hisat2.sam
else
   echo "Los datos ya estaban mapeados"
fi

## Pasar de .sam a .bam, ordenarlo e indexarlo
if [ ! -f ../results/mapeado/SRR13795616.sorted.bam ]
then
   samtools view -Sbh ../results/mapeado/SRR13795616_hisat2.sam > ../results/mapeado/SRR13795616_hisat2.bam
   samtools sort ../results/mapeado/SRR13795616_hisat2.bam -o ../results/mapeado/SRR13795616_hisat2.sorted.bam
   samtools index ../results/mapeado/SRR13795616_hisat2.sorted.bam
   ## Vemos unas estadísticas básicas
   samtools flagstat ../results/mapeado/SRR13795616_hisat2.sorted.bam > ../results/mapeado/SRR13795616_hisat2.sorted.flagstat
else
    echo "Los datos ya se habían pasado a bam y se habían ordenado e indexado"
fi
