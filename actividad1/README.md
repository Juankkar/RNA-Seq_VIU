# **Actividad 1 *From reads to gene counts***

## Objetivo: 

El objetivo de esta actividad es que el estudiante demuestre que ha adquirido las habilidades y competencias necesarias para llevar a cabo un correcto pre-procesamiento inicial de datos de RNA-seq, desde la caracterización y estudio de calidad de las lecturas, hasta su filtrado, mapeado y generación de la matriz de recuentos.

### Obtención de los datos:
La actividad consiste en el análisis de una muestra real de RNA-seq que forma parte de un estudio mayor sobre la eficacia del tratamiento con extracto de bacterias BV frente a la infección por betacoronavirus murino en un modelo de ratón.

### Directorios: 
* [code](code): código usado, conjunto de scripts y el ambiente usado para el análisis.
  * [enviroments](code/enviroments): ambiente usado para este análisis llamado [05MBIF](code/enviroments/05MBIF.yml).
  * [1descarga_raw](code/1descarga_raw.sh): descarga de los datos crudos a usar, los cuales serán depositados en la carpeta raw dentro de data.
  * [2_profundidad_len_fastqc](code/2_profundidad_len_fastqc.sh): en este script se realiza el cálculo manual de la profundidad de la secuenciación, y el rango de logitudes de las secuencias, para posteriormente realizar un análisis de la calidad más exhaustivo mediante fastQC.
  * [3_instalar_TrimGalore](code/3_instalar_TrimGalore.sh): Se instala Trim Galore y se dan las instrucciones para su uso.
  * [4_trimmeado](code/trimmeado.sh): arreglamos las secuencias con Trim Galore.
  * [5_mapeado](code/5_mapeado.sh): Mapeamos las secuencias con Hisat2 frente al genoma de referencia de *Mus musculus*
  * [6_pasar_sam_bam](code/6_pasar_sam_bam.sh): pasamos el archivo .sam obtenido del anterior script y lo pasamos a un archvo .bam.
  * [7_anotaciones_alineamiento](code/7_anotaciones_alineamiento.sh): En este script además de la descarga y anotación de los alineamientos mediante el archivo GTF de *Mus musculus*, terminamos con la formación de la matriz de conteos.
  * [8_all_code](code/8_all_code.sh): es un script cuya intención es correr todos los scripts al unísono. Tiene una serie de ***limitaciones*** como por ejemplo hay que tener Trim Galore instalado (script 3), con lo que está en desarrollo.
  * [9_reset](code/9_reset.sh): sirve para resetear las estructuras de los directorios data y results.
 
 
 * [data](data): conjunto de directorios con los datos, conforme se realicen los scripts se van rellenando.
   * [annotation](data/annotation): directorio con el que se almacenará principalmente archivo GTF de *Mus musculus*.
   * [processed](data/processed): datos crudos procesados.
   * [raw](data/raw): datos crudod de partida.
   * [reference_genome](data/reference_genome): donde se almacenará el genoma de referencia.
  
 * [results](results): resultados varios, htmls, graficas, pero además la matriz de recuentos final.
   * [fastqc_results](results/fastqc_results): html análisis de fastQC.
   * [mapeado](results/mapeado): Resultados del alineamiento con el genoma de referencia.
   * [tables](results/tables): En este directorio se almacenará la matriz de recuento final (No ha sido ignorada por .gitignore, con lo que puede visualizarse).

  
