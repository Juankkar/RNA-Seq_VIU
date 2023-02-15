#!/usr/bin/env bash

## Pasar de .sam a .bam, ordenarlo e indexarlo
if [ ! -f ../results/mapeado/SRR13795616_hisat2.bam ]
then
   samtools view -Sbh ../results/mapeado/SRR13795616_hisat2.sam > ../results/mapeado/SRR13795616_hisat2.bam
   samtools sort ../results/mapeado/SRR13795616_hisat2.bam -o ../results/mapeado/SRR13795616_hisat2.sorted.bam
   samtools index ../results/mapeado/SRR13795616_hisat2.sorted.bam
   ## Vemos unas estadísticas básicas
   samtools flagstat ../results/mapeado/SRR13795616_hisat2.sorted.bam > ../results/mapeado/SRR13795616_hisat2.sorted.flagstat
else
    echo "Los datos ya se habían pasado a bam y se habían ordenado e indexado"
fi

