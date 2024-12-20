options(stringsAsFactors = F)
library(tidyr)
library(clusterProfiler)
library(msigdbr)  #
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
library(ggplot2)
library(ggcor)
library(IOBR)
library(stringr)
library(ggExtra)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## Input your cancer type and gene

gene <- "PRTFDC1"
ca <- "TGCT"
### TCGA
## TCGA
load(file = "TCGA-TGCT-mrna_expr+pheno+sur.Rdata")
colnames(expr) <- gsub("\\.","-",colnames(expr))
## Annotation, change to Symbol
TGCT <- expr %>%  anno_eset(annotation = anno_grch38, probe = "id")
##
# cibersort <- deconvo_tme(eset = TGCT, method = "cibersort", arrays = FALSE, perm = 200)
# epic <- deconvo_tme(eset = TGCT, method = "epic", arrays = FALSE)
# mcp <- deconvo_tme(eset = TGCT, method = "mcpcounter")
# xcell <- deconvo_tme(eset = TGCT, method = "xcell", arrays = FALSE)
# timer <- deconvo_tme(eset = TGCT, method = "timer", group_list = rep(ca, dim(TGCT)[2]))
# estimate <- deconvo_tme(eset = TGCT, method = "estimate")
# quantiseq <- deconvo_tme(eset = TGCT, tumor = TRUE, arrays = TRUE, scale_mrna = TRUE, method = "quantiseq")
### Combine
# tme_combine <- cibersort %>% 
#   inner_join(., mcp, by = "ID") %>% 
#   inner_join(., xcell, by = "ID") %>%
#   inner_join(., epic, by = "ID") %>% 
#   inner_join(., timer, by = "ID") %>% 
#   inner_join(., estimate, by = "ID") %>% 
#   inner_join(., quantiseq, by = "ID")
# dim(tme_combine)
##
# save(tme_combine, file = paste0("TCGA-", ca, "_TME_combine.Rdata"))
load(file = paste0("TCGA-", ca, "_TME_combine.Rdata"))
TCGA_TME <- tme_combine %>% 
  column_to_rownames(var = "ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rbind(., TGCT[gene,]) %>% 
  .[grepl(paste0(("CD8|NK_|Tregs|Macrophage|Th1|Fibroblast|Dendritic|^DC_xCell|T_cells_CD8_CIBERSORT|CD8+_Tcm|CD8+_Tem|DC_TIMER|dendritic|NKcells|"), gene), rownames(.)),] %>% 
  .[!grepl("_naive_T-cells_xCell|Macrophages_xCell|_T-cells_xCell", rownames(.)),]

genes <- TGCT[gene,]

### Butterfly plot
TCGAimmSCORE.score <- TCGA_TME
TCGAimmCor <- NULL
for (i in rownames(TCGAimmSCORE.score)) {
  cr <- cor.test(as.numeric(TCGAimmSCORE.score[i,]),
                 as.numeric(genes),
                 method = "spearman")
  TCGAimmCor <- rbind.data.frame(TCGAimmCor,
                                 data.frame(gene = gene,
                                            path = i,
                                            r = cr$estimate,
                                            p = cr$p.value,
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
}
TCGAimmCor$spearman_cor <- round(TCGAimmCor$r, 2)
TCGAimmCor$method <-  str_match(TCGAimmCor$path, ".*_([a-zA-Z]*)")[,2]  %>% gsub("_", " ", .) 
TCGAimmCor$sign <- ifelse(TCGAimmCor$r > 0, "pos", "neg")
TCGAimmCor$absR <- abs(TCGAimmCor$r)
TCGAimmCor$rSeg <- as.character(cut(TCGAimmCor$absR, c(0, 0.25, 0.5, 0.75, 1), labels = c("0.25", "0.50", "0.75", "1.00"), include.lowest = T))
TCGAimmCor$pSeg <- as.character(cut(TCGAimmCor$p, c(0, 0.001, 0.01, 0.05, 1), labels = c("<0.001", "<0.01", "<0.05", "ns"), include.lowest = T))
TCGAimmCor[nrow(TCGAimmCor), "pSeg"] <- "Not Applicable"
TCGAimmCor$rSeg <- factor(TCGAimmCor$rSeg, levels = c("0.25", "0.50", "0.75", "1.00"))
TCGAimmCor$pSeg <- factor(TCGAimmCor$pSeg, levels = c("<0.001", "<0.01", "<0.05", "Not Applicable", "ns"))
TCGAimmCor$sign <- factor(TCGAimmCor$sign, levels = c("pos", "neg"))

for (i in 1:length(rownames(TCGAimmCor))) {
  if (TCGAimmCor$p[i] < 0.001) {
    TCGAimmCor$p_value[i] = "***"
  }
  if (TCGAimmCor$p[i] < 0.01 & TCGAimmCor$p[i] > 0.001) {
    TCGAimmCor$p_value[i] = "**"
  }
  if (TCGAimmCor$p[i] < 0.05 & TCGAimmCor$p[i] > 0.01) {
    TCGAimmCor$p_value[i] = "*"
  }
  if (TCGAimmCor$p[i] > 0.05) {
    TCGAimmCor$p_value[i] = round(TCGAimmCor$p[i], 2)
  }
}
library(xlsx)
write.xlsx(TCGAimmCor, file = paste(ca, gene, "immCor.xlsx", sep = "-"), row.names = F, sheetName = paste("TCGA", gene, sep = "-"))

# Single plot
### corplot
library(ggstatsplot)

sincorplot <- function(a, b, method = "spearman", title = "TCGA", subtitle = str_match(b, ".*_([a-zA-Z]*)")[,2]  %>% gsub("_", " ", .)) {
  
  plot_df <-  as.data.frame(t(TCGA_TME))[,c(a, b)]
  
  names(plot_df) <- c("A", "B")
  corT = cor.test(plot_df$A, plot_df$B, method = "spearman")
  cor = corT$estimate
  pValue = corT$p.value
  p1 = ggplot(plot_df, aes(A, B)) + 
    xlab(a) +
    ylab(b) +
    ggtitle(subtitle) +
    labs(x = paste0("Expression of ", gene, " "), y = "CD8+ T Cell infiltration level") +
    geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x) + 
    theme_classic() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          # plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold.italic")) +
    scale_y_continuous(labels = function(b) sprintf("%.2f", b)) +
    stat_cor(method = 'spearman', aes(x = A, y = B))
  p2 = ggMarginal(p1, type = "density", xparams = list(fill = "orange"), yparams = list(fill = "blue"))
  
  pdf(file = paste("TCGA ", b, ".pdf"), width = 4, height = 3)
  print(p2) 
  # export::graph2pdf(file = paste("TCGA ", b, ".pdf"), width = 3.5, height = 3.5)
  dev.off()
}

sincorplot(gene,"CD8+_Tcm_xCell")
sincorplot(gene,"CD8+_Tem_xCell")
sincorplot(gene,"CD8_T_cells_MCPcounter")
sincorplot(gene,"T_cell_CD8_TIMER")
sincorplot(gene,"CD8_Tcells_EPIC")
sincorplot(gene,"T_cells_CD8_CIBERSORT")












