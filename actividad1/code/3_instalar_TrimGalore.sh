#!/usr/bin/env bash

if [ ! -d  TrimGalore-0.6.10 ]
then
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
    tar xvzf trim_galore.tar.gz
    rm trim_galore.tar.gz
else
    echo "Trim Galore ya esta en el directorio"
fi

echo "Instrucciones:"
echo '1. Realizar en la carpeta creada: export PATH=$PATH:$PWD'
echo "2. trim_galore archivo.fastq.gz -o directorio_del_output"

