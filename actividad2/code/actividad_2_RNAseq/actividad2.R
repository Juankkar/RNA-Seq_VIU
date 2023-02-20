library(tidyverse)
library(edgeR)
library(cowplot)
library(glue)
library(ggtext)
library(venn)
library(pheatmap)
library(matrixStats)
library(ggrepel)
library(UpSetR)

#######################
## Pregunta número 1 ##
#######################
## Cargamos los metadatos ##
metadata <- read.delim("../../data/experimento_GSE167749.tsv", row.names = 1)
metadata

tabla <- data.frame(table(metadata$GEO, metadata$disease_state))
tabla %>% 
  mutate(Var2=factor(Var2,
                     levels=c("healthy", "infected", "infected+treated"),
                     labels=c("Sanos", "Infectados", "Infectados tratados"))) %>% 
  ggplot(aes(Var2,Var1, fill=Freq)) +
  geom_tile(show.legend = F) +
  labs(
    title = "Muestras frente a los grupos",
    y = "Muestras",
    x = "Grupos de tratameinto"
  ) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = .5, face="bold"))

grupos <- factor(metadata$disease_state)
table(grupos)

diseño <- model.matrix(~0+grupos)
colnames(diseño) <- c("healthy","infected","infected_treated")
rownames(diseño) <- rownames(metadata)

#######################
## Pregunta número 2 ##
#######################
recuentos <- read.delim("../../data/recuento_GSE167749.tsv",row.names = 1)
colnames(recuentos) <- row.names(metadata)

## Representación de los recuentos de la media frente la variaza
infected_counts <- recuentos %>% 
  mutate(mean_counts=apply(recuentos[,6:10],1,mean),
         variance_counts=apply(recuentos[,6:10],1,var)) %>%
  select(-starts_with(c("E", "C", "ET")))

infected_counts %>%
  ggplot(aes(mean_counts, variance_counts)) +
  geom_point(alpha=.5) +
  geom_abline(intercept = 0,slope = 1, color="red") +
  labs(
    title = "Relación de la media frente la varianza por gen",
    subtitle = "Muestra: Ratones infectados no tratados",
    y = "Recuento varianzas",
    x = "Recuento medias"
  ) +
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black")
  )

## Cuantos genes figuran en la matriz de recuentos? cauntos son por betacoronavirus?
# Un total de 55431
recuentos %>% nrow() 
# Un total de 55421 gene de ratón y 10 genes de beta coronavirus
recuentos %>%
  mutate(genes = rownames(recuentos),
         especie=case_when(grepl("ENSMUS", genes)~"Ratón",
                           !(grepl(("ENSMUS"), genes))~"Betacoronavirus")) %>%
  group_by(especie) %>%
  summarise(n=n())

# Realizamos los filtros:
filtrado1 <- recuentos[rowSums(recuentos) > 0,]
filtrado2 <- recuentos[rowSums(recuentos) >= 4,]
filtrado3 <- recuentos[rowSums(recuentos) >= 10,]

## Funcion para realizar el histograma mas eficiente
gghist <- function(df, vector, expresion){
  cuerpo=ggplot(data=df, aes(vector)) +
    geom_histogram(stat = "bin", 
                   bins = 25,
                   color="black",
                   fill="lightgray") + 
    labs(title = glue("Filtrado de expresión {expresion}"),
         x = "Media log(CPM)",
         y = "Número de genes") +
    scale_y_continuous(expand = expansion(0),
                       limits = c(0,15000)) +
    theme_test() +
    theme(
      plot.title = element_text(face = "bold", hjust = .5, size = 14),
      axis.title = element_text(face = "bold", size = 12))
  return(cuerpo)
  }
hist1 <- gghist(filtrado1, aveLogCPM(filtrado1), "> 0")
hist2 <- gghist(filtrado2, aveLogCPM(filtrado2), ">= 4")
hist3 <- gghist(filtrado3, aveLogCPM(filtrado3), ">= 10")
## Juntamos todos los histogramas
plot_grid(hist1,hist2,hist3, nrow = 1)

#######################
## Pregunta número 3 ##
#######################
preDGE <- filtrado2
# Creamos el objeto DGEList
y <- DGEList(preDGE,group = grupos)
y <- calcNormFactors(y,method = "TMM")

## Apartado 1 vemos los norm.factors < 1
y$samples %>%
  as_tibble() %>%
  filter(norm.factors < 1)

## Segundo apartado
no_normalizado <- cpm(y,normalized.lib.sizes = F, log=TRUE) %>%
  as.data.frame() %>% 
  mutate(genes = row.names(y$counts)) %>% 
  pivot_longer(-genes, names_to = "muestras", values_to = "valores") %>% # pasamos a formato tidy
  mutate(grupo=str_remove(pattern = "[:digit:]",muestras))

boxplot1 <- no_normalizado %>% 
  ggplot(aes(muestras, valores,fill=grupo)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = median(no_normalizado$valores), color="red") +
  labs(title = "Tamaños librerías (Sin normalizado)",
       x="muestras",
       y="Log~2~ CPM") +
  scale_fill_manual(values = c("white", "gray", "orange")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.y = element_markdown(),
    axis.text = element_text(color = "black")
  )

cpm <- cpm(y,normalized.lib.sizes = T, log=TRUE) 
df_cpm <- cpm %>%
  as.data.frame() %>% 
  mutate(genes = row.names(cpm)) %>% 
  pivot_longer(-genes, names_to = "muestras", values_to = "valores") %>% # pasamos a formato tidy
  mutate(grupo=str_remove(pattern = "[:digit:]",muestras))

boxplot2 <- df_cpm %>% 
  ggplot(aes(muestras, valores,fill=grupo)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = median(df_cpm$valores), color="red") +
  labs(title = "Tamaños librerías (Normalizado)",
       x="muestras",
       y="Log~2~ CPM") +
  scale_fill_manual(values = c("white", "gray", "orange")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.y = element_markdown(),
    axis.text = element_text(color = "black")
  )

plot_grid(boxplot1,boxplot2)

## Tercer apartado
plotMDS(y, col = as.numeric(y$samples$group), 
        gene.selection = "common")

#######################
## Pregunta número 4 ##
#######################
# Calcula los tres tipos de dispersión en los datos del RNA-seq problema 
# mediante la función estimateDisp y observa cómo la tagwise dispersion y 
# la trended dipersion
y <- estimateDisp(y, diseño, robust=TRUE)

## Creamos un tibble con los vectores que nos interes y filtramos los genes
## de interés
tibble(genes=rownames(y$counts), 
       tagwise_dispersion = y$tagwise.dispersion, 
       trended_dispersion = y$trended.dispersion) %>% 
  filter(genes %in% c("ENSMUSG00000078354", 
                      "ENSMUSG00000049176",
                      "55060280"))
## Segundo apartado
fit <- glmQLFit(y, diseño, robust = TRUE)
fit$coefficients %>% 
  as_tibble() %>% 
  mutate(genes=row.names(y$counts)) %>% 
  filter(genes %in% c("ENSMUSG00000078354", 
                      "ENSMUSG00000049176",
                      "55060280"))

#######################
## Pregunta número 5 ##
#######################
# Creamos los grupos
I_H_contrast <- makeContrasts(infected-healthy, levels=diseño)
IT_I_contrast <-makeContrasts(infected_treated-infected, levels=diseño)

IvsH <- glmQLFTest(fit, contrast = I_H_contrast)
ITvsI  <- glmQLFTest(fit, contrast = IT_I_contrast)

IvsH_DGE <- topTags(IvsH, n=Inf)
ITvsI_DGE <- topTags(ITvsI, n=Inf)

## Primer apartado ##
par(mfrow=c(2, 2)) ## Nos permitirá juntar los histogramas

hist(IvsH_DGE$table$PValue, 
     main = "Infectados vs Sanos", xlab = "p-Value")
hist(ITvsI_DGE$table$PValue, 
     main = "Infectados vs Sanos", xlab = "FDR")
hist(IvsH_DGE$table$FDR, 
     main = "Enfermos tratados vs Enfermos", xlab = "p-Value")
hist(ITvsI_DGE$table$FDR, 
     main = "Enfermos tratados vs Enfermos", xlab = "FDR")

## Segundo apartado ##
tibble(
Infected_vs_Healthy = IvsH_DGE$table %>% filter(FDR <= 0.05) %>% nrow(),
InfectesT_vs_Infected = ITvsI_DGE$table %>% filter(FDR <= 0.05) %>% nrow()
)

tibble(
  Infected_vs_Healthy = IvsH_DGE$table %>% filter(FDR <= 0.01) %>% nrow(),
  InfectesT_vs_Infected = ITvsI_DGE$table %>% filter(FDR <= 0.01) %>% nrow()
) 

## Hacemos el plot MDS
is.deIH <- decideTestsDGE(IvsH)
plotMD(IvsH, status=is.deIH, cex = 0.5, legend = "bottomright")
abline( h = c( -1, 1 ), col = "black")

is.deITI <- decideTestsDGE(ITvsI)
plotMD(ITvsI, status=is.deITI, cex = 0.5, legend = "bottomright")
abline( h = c( -1, 1 ), col = "black")

## Tercer apartado
IvsH_DGE_volcano <- IvsH_DGE$table %>% 
  mutate(genes_id=rownames(IvsH_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                             FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                             !((FDR <= 0.05 & logFC >=1)) & 
                               !((FDR <= 0.05 & logFC <=-1)) ~ "No UP/DOWN"))
IvsH_DGE_volcano %>% 
  group_by(expresion) %>% 
  summarise(n=n())

ITvsI_DGE_volcano <- ITvsI_DGE$table %>% 
  mutate(
    genes_id=rownames(ITvsI_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                              FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                              !((FDR <= 0.05 & logFC >=1)) & 
                                !((FDR <= 0.05 & logFC <=-1)) ~ "No UP/DOWN"))
ITvsI_DGE_volcano %>% 
  group_by(expresion) %>%
  summarise(n=n())
 

## Infectados vs Sanos
top_3_genes_IvsH <- IvsH_DGE$table %>% 
  arrange(PValue) %>% 
  head(n=3)

volcano1 <- IvsH_DGE_volcano %>% 
  mutate(expresion=factor(expresion,
                          levels=c("No UP/DOWN", "DOWN", "UP")),
         top_genes_labels = ifelse(genes_id %in% rownames(top_3_genes_IvsH),
                                   genes_id, NA)) %>% 
  ggplot(aes(logFC, -log10(FDR), color=expresion, label=top_genes_labels)) +
  geom_point(alpha=.5) +
  labs(
    title = "Volcano plot, DGE Ratones",
    subtitle = "Infectados vs sanos",
    x = "Log~2~ Fold Change",
    y = "-Log~10~ *P*",
    color="Expresión"
  ) +
  scale_color_manual(values = c("gray", "blue", "forestgreen")) +
  geom_vline(xintercept = c(1,-1), color="red", linetype="dashed") +
  geom_hline(yintercept = -log10(.05), color="red", linetype="dashed") +
  geom_text_repel(max.overlaps = Inf, show.legend = F, color="black") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = .5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

## Infectados tratados vs infectados
top_3_genes_ITvsI <- ITvsI_DGE$table %>% 
  arrange(FDR) %>% 
  head(n=3)

volcano2 <- ITvsI_DGE_volcano %>% 
  mutate(expresion=factor(expresion,
                          levels=c("No UP/DOWN", "DOWN", "UP")),
         top_genes_labels = ifelse(genes_id %in% rownames(top_3_genes_ITvsI),
                                   genes_id, NA)) %>% 
  ggplot(aes(logFC, -log10(FDR), color=expresion, label=top_genes_labels)) +
  geom_point(alpha=.5) +
  labs(
    title = "Volcano plot, DGE Ratones",
    subtitle = "Infectados tratados vs infectados",
    x = "Log~2~ Fold Change",
    y = "-Log~10~ *P*",
    color="Expresión"
  ) +
  scale_color_manual(values = c("gray", "blue", "forestgreen")) +
  geom_vline(xintercept = c(1,-1), color="red", linetype="dashed") +
  geom_hline(yintercept = -log10(.05), color="red", linetype="dashed") +
  geom_text_repel(max.overlaps = Inf, show.legend = F, color="black") +
  theme_classic() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = .5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown())

plot_grid(volcano1, volcano2)

## Último apartado
IvsH_DGE_up <- row.names(IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & IvsH_DGE$table$logFC>=1,])
IvsH_DGE_down <- row.names(IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & IvsH_DGE$table$logFC<=-1,])

ITvsI_DGE_up <- row.names(ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & ITvsI_DGE$table$logFC>=1,])
ITvsI_DGE_down <- row.names(ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & ITvsI_DGE$table$logFC<=-1,])

# Vennplot
venn(list("IvsH_up"=IvsH_DGE_up, "IvsH_down"=IvsH_DGE_down, 
          "ITvsI_up"=ITvsI_DGE_up, "ITvsI_down"=ITvsI_DGE_down ),
           zcolor = c("#999999", "#E69F00", "#56B4E9", "#009E73"))


## UpSet plot
lista_genes <- list(
  IvsH_DGE_up = IvsH_DGE_up,
  IvsH_DGE_down = IvsH_DGE_down,
  ITvsI_DGE_up = ITvsI_DGE_up,
  ITvsI_DGE_down = ITvsI_DGE_down
)

upset(fromList(lista_genes), order.by = "freq")

#######################
## Pregunta número 6 ##
#######################
## Primer apartado:
IvsH_DGE_filt <- IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & abs(IvsH_DGE$table$logFC)>=1,]
ITvsI_DGE_filt <- ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & abs(ITvsI_DGE$table$logFC)>=1,]
all_DGE_filt <- intersect(row.names(IvsH_DGE_filt), row.names(ITvsI_DGE_filt))
# Extracción de datos normalizados
cpms <- cpm(y$counts, log =TRUE)
colnames(cpms) <- rownames(metadata)
head(cpms)
# Hacemos un subset con los genes de interes
matrix_DGE_cpm <- cpms[all_DGE_filt,]
head(rowSds(matrix_DGE_cpm))
# Transformar la matriz de datos a Z-Score (función scale)
zscore_fun <- function(x){(x-mean(x))/sd(x)}
# Y escalamos
matrix_DGE_zscore <- t(apply(matrix_DGE_cpm,1,zscore_fun))
head(rowSds(matrix_DGE_zscore))

matrix_DGE_cpm_long <- matrix_DGE_cpm %>%
  as.data.frame() %>%
  mutate(matriz = rep("Matriz CPM", nrow(matrix_DGE_cpm))) %>% 
  pivot_longer(-matriz, names_to = "muestras", values_to = "valores")

matrix_DGE_zscore_long <- matrix_DGE_zscore %>%
  as.data.frame() %>%
  mutate(matriz = rep("Matriz Z-score", nrow(matrix_DGE_cpm))) %>% 
  pivot_longer(-matriz, names_to = "muestras", values_to = "valores")

rbind(matrix_DGE_cpm_long, matrix_DGE_zscore_long) %>%
  ggplot(aes(valores, fill=matriz)) +
  geom_density(color="black", alpha=.75) +
  labs(title = "Antes y después de realizar Z-score",
       x = "Valores") +
  scale_x_continuous(limits = c(-4, 12),
                     breaks = seq(-4,12,2)) +
  scale_fill_manual(name=NULL,
                    values = c("orange", "black")) +
  theme_test()+ 
  theme(legend.position = c(.8,.7))

plot(density(matrix_DGE_cpm))
plot(density(matrix_DGE_zscore))


for(i in c("euclidean", "manhattan")){
  pheatmap(
    matrix_DGE_cpm,
    cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
    cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
    show_rownames = FALSE, # There are too many genes to clearly show the labels,
    clustering_distance_rows = i,
    clustering_distance_cols = i,
    main = glue("Heatmap intersection DEG {i}"),
    colorRampPalette(c(
    "skyblue",
              "white",
              "red"
    ))(25
    ),
    scale = "row") # Scale values in the direction of genes (rows)
}


## Último apartado, filtramos el coronavirus
madtrix_betacoronavirus <- matrix_DGE_cpm %>%
  as.data.frame() %>% 
  filter(str_detect(rownames(matrix_DGE_zscore), "^\\d")) %>%
  as.matrix()

for(i in c("euclidean", "manhattan")){
pheatmap(
  madtrix_betacoronavirus,
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_rownames = FALSE, 
  clustering_distance_rows = i,
  clustering_distance_cols = i,
  main = glue("Heatmap intersection DEG {i}"),
  method="complete",
  colorRampPalette(c(
    "skyblue",
             "white",
             "red"
  ))(25
  ),
  scale = "row")
}
#######################
## Pregunta número 7 ##
#######################
library(GO.db) ## OJO, enmascara funciones de dplyr, no bueno asi que lo dejo por aqui
require(org.Mm.eg.db)

library(limma)
base <- as.data.frame(org.Mm.egENSEMBL2EG)
IvsH_DGE_up <- row.names(IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & IvsH_DGE$table$logFC>=1,])
IvsH_DGE_down <- row.names(IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & IvsH_DGE$table$logFC<=-1,])
ITvsI_DGE_up <- row.names(ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & ITvsI_DGE$table$logFC>=1,])
ITvsI_DGE_down <- row.names(ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & ITvsI_DGE$table$logFC<=-1,])

buscar_go <- function(names_regulated){
  cuerpo = for(i in c("BP", "MF")) {
    print(glue("===> {i} <==="))
    funer <- base[base$ensembl_id %in% names_regulated, ]
    go <- goana(funer$gene_id, species="Mm", geneid= "ENTREZID") # Busqueda de go y estudio de su frecuencia
    print(topGO(go, ontology=i, number = 3))
  }
  return(cuerpo)
}

buscar_go(IvsH_DGE_up)
buscar_go(IvsH_DGE_down)
buscar_go(ITvsI_DGE_up)
buscar_go(ITvsI_DGE_down)

## Segundo apartado ##
## Uso de DOSE
library(clusterProfiler)
funer <- base[base$ensembl_id %in% IvsH_DGE_up, ]
kegg <- enrichKEGG(gene = funer$gene_id,
                       organism = "mmu",
                       pvalueCutoff = .05,
                       qvalueCutoff = .05)
barplot(kegg)
