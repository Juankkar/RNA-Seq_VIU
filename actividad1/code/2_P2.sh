#!/usr/bin/env bash

echo ""
echo "############################################"
echo "EJERCICIO1: PROFUNDIDAD, LONGITUD SECUENCIAS"
echo "############################################"
echo ""

echo ">>> Número de secuencias del formato .fastq"
num_linea=$(echo $(awk {'print NR'} ../data/raw/SRR13795616.fastq.gz | tail -1))
echo $((num_linea/4))

## Para el promedio de longitud Filtramos las secuencias del archivo
zcat ../data/raw/SRR* | sed -n '2~4p' > ../data/processed/ejercicio1/secuencias.txt
## Para ello usaremos un bucle while y haremos un pipe a un awk para realizar el 
## cálculo de la media de longitud de las secuencias
echo ">>> Longitud promedio del las secuencias: "
while read -r line
do 
    echo ${#line}
done < ../data/processed/ejercicio1/secuencias.txt | awk 'BEGIN{avg=0}{avg += $1}END{print avg/NR}' 

echo ">>> Finalmente realizamos un fastqc, en caso de no haberlo hecho"
if [ ! -f ../results/fastqc_results/SRR13795616_fastqc.html ]
then
    fastqc ../data/raw/SRR13795616.fastq.gz
    mv ../data/raw/*fastqc* ../results/fastqc_results
else
    echo "Los datos ya se encuentran en ../results/fastqc_results"
fi

