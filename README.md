# Actividades para realizar RNA-seq Universidad Internacional de Valencia (VIU)

## 2 actividades originales:

* [Actividad 1](actividad1): ***From reads to gene counts***: preprocesado de datos de RNA-seq para la obtención de una matriz de recuentos.
* [Actividad 2](actividad2): ***From genecounts to DEG and pathways***. Análisis estadísticos para genes expresados diferencialmente (DGE).

Ahora las he convertido en un proyecto.


## Directorios

### [Code](code): directorio con los ambientes y programas desarrollados.

* Tenemos por un lado otro directorio con los ambientes a usar: [environments](code/environments/).

    * [data_analysisR](code/environments/data_analysisR.yml): este ambiente sirve para realizar estudios estadísticos con el lenguaje de programación R, para el R markdown de la última regla del archivo Snakefile.

    * [htseq](code/environments/htseq.yml): software necesario para realizar el conteo de los genes expresados, tras el mapeo.

    * [rnaseq](code/environments/rnaseq.yml): mayoría del software para realizar RNA-Seq, que no requieran distintas dependencias.

    * [trimming](code/environments/trimming.yml): ambiente necesario para usar TirmGalore!. Se encuentra por separado, por temas de dependencias.

* Por otro tenemos los programas desarrollados.

    * [1download_data.sh](code/1download_data.sh): descarga los datos de un repositrio público en el ENA.
    * [2num_lineas.py](code/2num_lineas.py): script de python para realziar una de las actividades pediadas por la profesora.
    * [3anotaciones.sh](code/3anotaciones.sh): descargamos la base de datos de las anotaciones del genoma de la especie *Mus musculus*.
    * [4data_analysis.Rmd](code/4data_analysis.Rmd): R markdown que permite realizar un informe básico para resolver lo que sería la actividad 2, que consiste en el análisis de datos transcriptómicos obtenidos.

### [Data](data/): directorio con los datos de estidio.

* [annotations](data/annotations/): donde se almacena la base de datos de las anotaciones del genoma de *Mus musculus*.

* [processed](data/processed/): Secuencias procesadas con TinGalore!.

* [raw](data/raw/): Secuencias crudas salidas del secuenciador.

* [reference](data/reference/): genoma de referencia de *Mus musculus* 

### [Metadata](metadata/): directorio con alguno de los metadatos u otros porporcionados por la profesora.

* [betacoronavirus.tsv](metadata/betacoronavirus.tsv): conteo de los genes de Betacoronavirus proporcionados por la profesora.

* [experimento_GSE167749.tsv](metadata/experimento_GSE167749.tsv): metadatos del proyecto original, proporcionados por la profesora.

* [report.tsv](metadata/report.tsv): metadatos del proyecto, con las muestras y los links para descargar los archivos, que serán usados en el script ***1download_data.sh***.

### [results](results/): directorio donde se irán almacenado los resultados.

### Archivos de Snakemake:

* [config.yaml](config.yaml): archivo que permite automatizar/dar más flexibilidad al archivo de Snakfile, con las muestras de estudio y parámetros.

* [Snakefile](Snakefile): Archivo con las reglas a seguir de forma jerárquica (solo será necesario ejecutarlas 1 vez cada una de ellas).

### Otros archivos

* [actividades_originales.zip](actividades_originales.zip): las actividades originales que fueron entregadas.

