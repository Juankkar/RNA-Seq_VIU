#!/usr/bin/env bash
echo ""
echo "############################################"
echo "EJERCICIO2: PROFUNDIDAD, LONGITUD SECUENCIAS"
echo "############################################"
echo ""

num_linea=$(echo $(awk {'print NR'} ../data/raw/SRR13795616.fastq.gz | tail -1))

if [ ! -f ../data/processed/ejercicio2/len_seqs.txt ]
then
    ## Para el rango de la longitud Filtramos las secuencias del archivo
    zcat ../data/raw/SRR* | sed -n '2~4p' > ../data/processed/ejercicio2/secuencias.txt
    ## Para ello usaremos un bucle while y calculando la longitud de cada secuencia
    while read -r line
        do 
        echo ${#line}
        done < ../data/processed/ejercicio2/secuencias.txt > ../data/processed/ejercicio2/len_seqs.txt 
else
    echo "Ya tenemos la longitud de las secuencias cargadas"
fi

max=$(sort ../data/processed/ejercicio2/len_seqs.txt | tail -n 1)
min=$(sort ../data/processed/ejercicio2/len_seqs.txt | head -n 1)

echo ""
echo "=-----------------  RESULTADOS    -------------------="
echo "====>       PROFUNDIDAD = $((num_linea/4))       <===="
echo "====>  VALOR MINIMO = $min | VALOR MAXIMO = $max <===="
echo "=----------------------------------------------------="
echo ""

echo ">>> Finalmente realizamos un fastqc, en caso de no haberlo hecho"
if [ ! -f ../results/fastqc_results/SRR13795616_fastqc.html ]
then
    fastqc ../data/raw/SRR13795616.fastq.gz
    mv ../data/raw/*fastqc* ../results/fastqc_results
else
    echo "Los datos ya se encuentran en ../results/fastqc_results"
fi

