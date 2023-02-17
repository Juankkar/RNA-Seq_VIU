library(tidyverse)
library(edgeR)
library(cowplot)
library(glue)
library(ggtext)
library(venn)
library(pheatmap)
library(matrixStats)

metadata <- read.delim("../../data/experimento_GSE167749.tsv", row.names = 1)
metadata

#######################
## Pregunta número 1 ##
#######################
heatmap(table(metadata$disease_state, metadata$GEO))

grupos <- factor(metadata$disease_state)
table(grupos)

diseño <- model.matrix(~0+grupos)
colnames(diseño) <- c("healthy","infected","infected_treated")
rownames(diseño) <- rownames(metadata)

recuentos <- read.delim("../../data/recuento_GSE167749.tsv",row.names = 1)
colnames(recuentos) <- row.names(metadata)

#######################
## Pregunta número 2 ##
#######################
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
    caption = "Lo que vemos es que cuando comparamos el valor 
    medio frente la varianza la relación no so iguales, esto se debe
    a que muchos genes se encuentran expresándose más que otros, sobre 
    todo, a valores promedio más altos de expresión más altos, existie
    más variacón den los datos.",
    y = "Recuento varianzas",
    x = "Recuenteo medias"
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
  as_tibble() %>%
  mutate(genes = rownames(recuentos)) %>%
  select(genes) %>%
  mutate(especie=case_when(grepl("ENSMUS", genes)~"Ratón",
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
## Primer apartado
# eliminamos la columna de las sumas que hicimos en un principio
preDGE <- filtrado2
# Creamos el objeto DGEList
y <- DGEList(preDGE,group = grupos)
y <- calcNormFactors(y,method = "TMM")


## Apartado 1 vemos los norm.factors < 1
y$samples %>%
  as_tibble() %>%
  filter(norm.factors < 1)

## Segundo apartado
no_normalizado <- y$counts %>%
  as.data.frame() %>% 
  mutate(genes = row.names(y$counts)) %>% 
  pivot_longer(-genes, names_to = "muestras", values_to = "valores") %>% # pasamos a formato tidy
  mutate(grupo=str_remove(pattern = "[:digit:]",muestras))

boxplot1 <- no_normalizado %>% 
  ggplot(aes(muestras, valores/10000,fill=grupo)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(yintercept = median(no_normalizado$valores)/10^4, color="red") +
  labs(title = "Tamaños librerías (Sin normalizar)",
       x="muestras",
       y="CPM x 10⁴") +
  scale_fill_manual(values = c("white", "gray", "orange")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
    axis.title = element_text(face = "bold", size = 12),
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
  labs(title = "Tamaños librerías (logCPMs)",
       x="muestras",
       y="Log2 CPM") +
  scale_fill_manual(values = c("white", "gray", "orange")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = .5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = .5, face = "italic", size = 12),
    axis.title = element_text(face = "bold", size = 12),
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
y$common.dispersion
head(y$trended.dispersion)
head(y$tagwise.dispersion)
head(y$AveLogCPM)

# Representa en un gráfico de la variabilidad biológica que existe en el conjunto de datos con la función plotBCV.
df_dispersiones <- tibble(genes=row.names(y$counts), 
       tawise_dispersion = y$tagwise.dispersion, 
       trended_dispersion = y$trended.dispersion) 

df_dispersiones %>% filter(genes %in% c("ENSMUSG00000078354", 
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

## Combinamos ambas comparaciones en un unico df tidy
nrow_IvsH_DGE_table <- nrow(IvsH_DGE$table)
IvsH_DGE_table <- IvsH_DGE$table %>% 
  mutate(grupos=rep("Infectados vs Sanos"), nrow(nrow_IvsH_DGE_table))

nrow_ITvsI_DGE_table <- nrow(ITvsI_DGE$table)
ITvsI_DGE_table <- ITvsI_DGE$table %>% 
  mutate(grupos=rep("Infectados tratados vs Infectados"), nrow(nrow_IvsH_DGE_table))

comparaciones_tidy <- rbind(IvsH_DGE_table, ITvsI_DGE_table) 

hist_DGE_pvalue <- comparaciones_tidy %>%
  as_tibble() %>%
  ggplot(aes(PValue, fill = grupos)) +
  geom_histogram(color="black", show.legend = FALSE) +
  geom_vline(xintercept = 0.05, color="red", linetype="dashed") +
  facet_grid(~grupos) +
  labs(title = "Histogramas de los p-Values",
       y = "Recuentos",
       x = "p-Value")

hist_DGE_fdr <- comparaciones_tidy %>%
  as_tibble() %>%
  ggplot(aes(FDR, fill = grupos)) +
  geom_histogram(color="black", show.legend = FALSE) +
  geom_vline(xintercept = 0.05, color="red", linetype="dashed") +
  facet_grid(~grupos) +
  labs(title = "Histogramas de los FDR",
       y = "Recuentos",
       x = "FDR")
plot_grid(hist_DGE_pvalue, hist_DGE_fdr, ncol = 1)

comparaciones_tidy %>%
  filter(grupos == "Infectados vs Sanos" & FDR >= 0.05 ) %>%
  nrow()

comparaciones_tidy %>% 
  filter(grupos == "Infectados vs Sanos" & FDR <= 0.01 ) %>% 
  nrow()

comparaciones_tidy %>% 
  filter(grupos == "Infectados tratados vs Infectados" & FDR <= 0.05 ) %>% 
  nrow()

comparaciones_tidy %>% 
  filter(grupos == "Infectados tratados vs Infectados" & FDR <= 0.01) %>% 
  nrow()

# Forma 2: Seleccionamons los genes con un FDR menor o igual de 0.05.
is.deIH <- decideTestsDGE(IvsH)
plotMD(IvsH, status=is.deIH, cex = 0.5, legend = "bottomright")
abline( h = c( -1, 1 ), col = "black")

is.deITI <- decideTestsDGE(ITvsI)
plotMD(ITvsI, status=is.deITI, cex = 0.5, legend = "bottomright")
abline( h = c( -1, 1 ), col = "black")

## Tercer apartado
IvsH_DGE$table %>% 
  mutate(genes_id=rownames(IvsH_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                             FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                             !((FDR <= 0.05 & logFC >=1)) & 
                               !((FDR <= 0.05 & logFC <=-1)) ~ "NS")) %>% 
  group_by(expresion) %>% 
  summarise(n=n())

ITvsI_DGE$table %>% 
  mutate(genes_id=rownames(IvsH_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                             FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                             !((FDR <= 0.05 & logFC >=1)) & 
                               !((FDR <= 0.05 & logFC <=-1)) ~ "NS")) %>% 
  group_by(expresion) %>%
  summarise(n=n())
 

## Infectados vs Sanos
top_3_genes_IvsH <- IvsH_DGE$table %>% 
  arrange(PValue) %>% 
  head(n=3)

volcano1 <- IvsH_DGE$table %>% 
  mutate(genes_id=rownames(IvsH_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                             FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                             !((FDR <= 0.05 & logFC >=1)) & 
                                 !((FDR <= 0.05 & logFC <=-1)) ~ "NS"),
         expresion=factor(expresion,
                          levels=c("NS", "DOWN", "UP")),
         top_genes_labels = ifelse(genes_id %in% rownames(top_3_genes_IvsH),
                                   genes_id, NA)) %>% 
  ggplot(aes(logFC, -log10(FDR), color=expresion, label=top_genes_labels)) +
  geom_point(alpha=.5) +
  labs(
    title = "Volcano plot, DGE Ratones",
    subtitle = "Infectados vs sanos",
    x = "Log~2~ Fold Change",
    y = "-Log~10~ *P*"
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
  arrange(PValue) %>% 
  head(n=3)

volcano2 <- ITvsI_DGE$table %>% 
  mutate(genes_id=rownames(ITvsI_DGE$table),
         expresion=case_when(FDR <= 0.05 & logFC >=1 ~ "UP",
                             FDR <= 0.05 & logFC <=-1 ~ "DOWN",
                             !((FDR <= 0.05 & logFC >=1)) & 
                               !((FDR <= 0.05 & logFC <=-1)) ~ "NS"),
         expresion=factor(expresion,
                          levels=c("NS", "DOWN", "UP")),
         top_genes_labels = ifelse(genes_id %in% rownames(top_3_genes_ITvsI),
                                   genes_id, NA)) %>% 
  ggplot(aes(logFC, -log10(FDR), color=expresion, label=top_genes_labels)) +
  geom_point(alpha=.5) +
  labs(
    title = "Volcano plot, DGE Ratones",
    subtitle = "Infectados tratados vs infectados",
    x = "Log~2~ Fold Change",
    y = "-Log~10~ *P*"
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

venn(list("IvsH_up"=IvsH_DGE_up, "IvsH_down"=IvsH_DGE_down, 
          "ITvsI_up"=ITvsI_DGE_up, "ITvsI_down"=ITvsI_DGE_down ),
           zcolor = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

#######################
## Pregunta número 6 ##
#######################
## Primer apartado:
IvsH_DGE_filt <- IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & abs(IvsH_DGE$table$logFC)>=1,]
ITvsI_DGE_filt <- ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & abs(ITvsI_DGE$table$logFC)>=1,]
all_DGE_filt <- intersect(row.names(IvsH_DGE_filt),row.names(ITvsI_DGE_filt))

cpms <- cpm(y$counts, log =TRUE)
colnames(cpms) <- rownames(metadata)

matrix_DGE_cpm <- cpms[all_DGE_filt,]
head(rowSds(matrix_DGE_cpm))

## Escalamos
# creamos nuestra propia función
zscore_fun <- function(x){(x-mean(x))/sd(x)}
# Y escalamos
matrix_DGE_zscore <- t(apply(matrix_DGE_cpm,1,zscore_fun))
head(rowSds(matrix_DGE_zscore))
