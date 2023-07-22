configfile: "config.yaml"


rule all:
  input:
    expand("results/visualizar_datos/{sample}.txt", 
           sample=config["samples"]),
    expand("results/fastqc/{sample}.html", 
               sample=config["samples"]),
    expand("data/processed/{sample}_trimming.fq.gz",
           sample=config["samples"]),
    expand("results/fastqc/trimmed/{sample}_trimmed_fastqc.html",
           sample=config["samples"]),
    expand("results/mapped_reads/{sample}_hisat2.sam",
           sample=config["samples"]),
    expand("results/mapped_reads/bam_files/{sample}_hisat2.bam",
           sample=config["samples"]),
    expand("results/mapped_reads/bam_files/{sample}_hisat2_sorted.bam",
           sample=config["samples"]),
    expand("results/visualizar_datos/pregunta7/{sample}.txt",
           sample=config["samples"]),
    expand("results/tables/{sample}_counts.tsv",
            sample=config["samples"]) 
 
#------------------------------------#
#             Set-up                 #
#------------------------------------#


rule download_data:
    input:
        script = "code/1download_data.sh"
    output:
        touch("task/download_data.done")
    params:
        field=config["field_download"]
    conda:
      "code/environments/rnaseq.yml"
    shell:
        """
        bash {input.script} {params.field}
        """


rule reference_genome:
    output: 
        "data/reference/genome.fa"
    params:
        genome=config["reference_genome"]
    conda:
      "code/environments/rnaseq.yml" 
    shell:
        """
        wget -O {output}.gz {params.genome} 
        gzip -d {output}.gz
        """


#-------------------------------------------#
#         Start of the Workflow             #
#-------------------------------------------#


rule ver_datos:
  input:
    data = "data/raw/{sample}.fastq.gz",
    script = "code/2num_lineas.py"
  output:
    "results/visualizar_datos/{sample}.txt"
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    """
    num_linea=$(echo $(zcat {input.data} | awk {{'print NR'}} | tail -1))

    zcat {input.data} | sed -n '2~4p' > results/visualizar_datos/temporal.txt

    python {input.script} > {output}

    max=$(tail -n 1 {output})
    min=$(head -n 1 {output})

    echo ""
    echo "=-----------------  RESULTADOS    -------------------="
    echo "====>       PROFUNDIDAD = $((num_linea/4))       <===="
    echo "====>  VALOR MINIMO = $min | VALOR MAXIMO = $max <===="
    echo "=----------------------------------------------------="
    echo ""

    rm results/visualizar_datos/temporal.txt
    """


rule fastqc:
  input:
    "data/raw/{sample}.fastq.gz"
  output:
    "results/fastqc/{sample}_fastqc.html"
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    "fastqc {input} -o results/fastqc"


rule trimming:
  input:
    "data/raw/{sample}.fastq.gz"
  output:
    "data/processed/{sample}_trimmed.fq.gz"
  conda:
    "code/environments/trimming.yml"
  shell:
    """
    trim_galore {input} -o data/processed/
    """

rule fastqc_trimmed:
  input:
    "data/processed/{sample}_trimmed.fq.gz"
  output:
    "results/fastqc/trimmed/{sample}_trimmed_fastqc.html"
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    "fastqc {input} -o results/fastqc/trimmed"


rule mapping:
  input:
    samples_trimmed = "data/processed/{sample}_trimmed.fq.gz",
    reference = "data/reference/genome.fa"
  output:
    "results/mapped_reads/{sample}_hisat2.sam"
  params:
    strand = config["strandness"]
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    """
    ## Indexamos el genoma
    hisat2-build {input.reference} data/reference/genome
    
    ## Hacemos el mapeado es importante que lo hagamos 
    #  identifiquemos aquellas muestras que son reverse 
    #  y forward 
    hisat2 -k1 -U {input.samples_trimmed} \
        -x data/reference/genome \
        --rna-strandness {params.strand} \
        -S {output}
    """


rule sam_to_bam:
  input:
    "results/mapped_reads/{sample}_hisat2.sam"
  output:
    "results/mapped_reads/bam_files/{sample}_hisat2.bam",
  conda:
    "code/environments/rnaseq.yml"
  shell:
    """
    ## Convertimos los archivos SAM a BAM 
    samtools view -Sbh {input} > {output} 
    """


rule sorting_bam:
  input:
    "results/mapped_reads/bam_files/{sample}_hisat2.bam"
  output:
    "results/mapped_reads/bam_files/{sample}_hisat2_sorted.bam",
  conda:
    "code/environments/rnaseq.yml"
  shell:
    """
    ## Ordeamos los archivos
    samtools sort {input} > {output}

    ## Indexamos los archivos
    samtools index {output} 

    ## Creamos unas estadísticas del mapeo para tener contexto
    samtools flagstats {output} > {output}.flagstats
    """

rule download_annotations:
  input:
    "code/3anotaciones.sh" 
  output:
    "data/annotations/genome.gtf"
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    """
    bash {input}
    """

rule pregunta7:
  output:
    "results/visualizar_datos/pregunta7/{sample}.txt"
  conda:
    "code/environments/rnaseq.yml" 
  shell:
    """
    ## Responder a las siguientes cuestiones:
    echo "===> Número de Bases de datos en SOURCE: <===" \
      > {output} 
    grep -v ^# data/annotations/genome.gtf \
      | cut -f2 \
      | sort \
      | uniq -c \
      | sort \
      >> {output}
    
    echo "====> Features totales <====" \
      >> {output}
    grep -v ^# data/annotations/genome.gtf \
      | cut -f3 \
      | sort \
      | uniq -c \
      | sort \
      >> {output}
    """


rule counts:
  input:
    annotations = "data/annotations/genome.gtf",
    bam = "results/mapped_reads/bam_files/{sample}_hisat2_sorted.bam"
  output:
    "results/tables/{sample}_counts.tsv"
  conda:
    "code/environments/htseq.yml"
  shell:
    """
    htseq-count -f bam \
      -r pos \
      -m union \
      -s reverse \
      -t exon \
      -i gene_id \
      {input.bam} {input.annotations} \
      > {output}
    """

rule data_analysisR:
  input:
    "code/4data_analysis.Rmd"
  output:
    "results/4data_analysis.html"
  params:
    dir1 = "results/data_analysis/",
    dir2 = "results/data_analysis/plots/",
    dir3 = "results/data_analysis/tables/"
  conda:
    "code/environments/data_analysisR.yml"
  shell:
    """
    if [[ ! -d results/data_analysis/ ]]
    then
      for dir in {params.dir1} {params.dir2} {params.dir3} 
      do
        mkdir $dir 
      done
    fi

    R -e "library(rmarkdown); render('{input}')"
    mv code/4data_analysis.html results/ 
    """