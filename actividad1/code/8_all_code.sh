#!/usr/bin/env bash

array=(1descarga_raw.sh 2_profundidad_len_fastqc.sh 4_trimmeado.sh 5_mapeado.sh 6_anotaciones_alienamientos.sh) 

for i in ${array[*]} 
do 
    bash $i 
done
