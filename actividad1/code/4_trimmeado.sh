#!/usr/bin/env bash

## Realizamos el procesado de los datos crudos
trim_galore ../data/raw/SRR13795616.fastq.gz -o ../data/processed/ejercicio2_trimmed

## Realizamos el procesado de los datos trimmeados y lo enviamos a la
## carpeta de los resultados
fastqc ../data/processed/ejercicio2_trimmed/*.fq.gz -o ../results/fastqc_results

