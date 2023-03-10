#!/usr/bin/env bash

#####################################################################
## Descargamos el archivo de las notaciones en caso de no ternerlo ##
#####################################################################
if [ ! -f ../data/annotation/genome.gtf ]
then 
    wget -P ../data/annotation https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
    gunzip ../data/annotation/Mus_musculus.GRCm38.102.gtf.gz
    mv ../data/annotation/Mus_musculus.GRCm38.102.gtf ../data/annotation/genome.gtf
    echo "Se han descargado los datos correctamente"
else
    echo "Los datos ya estaban descargados"
fi

## Respondemos a las cuestiones: 
## ¿Cuántos programas y bases de datos se emplearon para anotar este archivo? (source)
## ¿Cuáles y cuántas FEATURES podemos encontrar en el archivo de anotaciones?
echo "====> Número de Bases de datos en SOURCE: <====" > ../data/annotation/pregunta7.txt
grep -v ^# ../data/annotation/genome.gtf | cut -f2 | sort | uniq -c | sort >> ../data/annotation/pregunta7.txt
echo "====> Features totales <====" >> ../data/annotation/pregunta7.txt
grep -v ^# ../data/annotation/genome.gtf | cut -f3 | sort | uniq -c | sort >> ../data/annotation/pregunta7.txt



## Y en caso de no tener la matriz de recuentos ya creada lo hacemos
if [ ! -f ../results/tables/SRR13795616_counts.tsv ] && [ ! -f ../data/results/tables/counts_no_s.tsv ]
then
    htseq-count -f bam -r pos -m union -s reverse -t exon -i gene_id ../results/mapeado/SRR13795616_hisat2.sorted.bam ../data/annotation/genome.gtf > ../results/tables/SRR13795616_counts.tsv
    htseq-count -f bam -r pos -m union -t exon -i gene_id ../results/mapeado/SRR13795616_hisat2.sorted.bam ../data/annotation/genome.gtf > ../results/tables/counts_no_s.tsv
    echo "Se ha creado una matriz de recuentos"
else
    echo "La matriz de recuentos ya estaba creada"
fi

