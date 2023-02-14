#!/usr/bin/env bash

array=(1descarga_raw.sh 2_P2.sh 3_instalar_TrimGalore.sh 4_trimmeado.sh 5_mapeado.sh 6_anotaciones_alienamientos.sh)

for i in ${array[*]}
do
    bash $i
done
