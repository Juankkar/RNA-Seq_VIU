metadata <- read.delim("Result/Tables/experimento_GSE167749.tsv",row.names = 1)
grupos <- factor(metadata$disease_state)
table(grupos)

diseño <- model.matrix(~0+grupos)
colnames(diseño) <- c("healthy","infected","infected_treated")
rownames(diseño) <- rownames(metadata)

recuentos <- read.delim("Result/Tables/recuento_GSE167749.tsv", row.names = 1)
recuentos <- recuentos[rowSums(recuentos) > 0, ]

y <- DGEList(recuentos,group = grupos)
y <- calcNormFactors(y,method = "TMM")
y <- estimateDisp(y, diseño, robust=TRUE)
fit <- glmQLFit(y, diseño, robust=TRUE)

I_H_contrast <- makeContrasts(infected-healthy, levels=diseño)
IT_I_contrast <-makeContrasts(infected_treated-infected, levels=diseño)

IvsH <- glmQLFTest(fit, contrast = I_H_contrast)
ITvsI  <- glmQLFTest(fit, contrast = IT_I_contrast)

IvsH_DGE <- topTags(IvsH, n=Inf)
ITvsI_DGE <- topTags(ITvsI, n=Inf)

write.table(IvsH_DGE$table, "Result/Tables/Infected_vs_Control_DGE.tsv", quote = F)
write.table(ITvsI_DGE$table, "Result/Tables/Infected_Treated_vs_Infected_DGE.tsv",quote = F)

IvsH_DGE_filt <- IvsH_DGE[IvsH_DGE$table$FDR<=0.05 & abs(IvsH_DGE$table$logFC)>=1,]
nrow(IvsH_DGE_filt)

ITvsI_DGE_filt <- ITvsI_DGE[ITvsI_DGE$table$FDR<=0.05 & abs(ITvsI_DGE$table$logFC)>=1,]
nrow(ITvsI_DGE_filt)

IvsH_DGE$table %>% 
  as_tibble() %>%
  mutate(expresion=case_when())
  ggplot(aes(logFC, -log10(FDR))) +
  geom_point()
