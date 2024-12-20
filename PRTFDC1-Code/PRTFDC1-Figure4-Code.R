Here is the translated version of your R script with comments and text in English:
  
  ```r
rm(list = ls())
## Load R packages
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(Biobase)
library(BiocGenerics)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(limma) # Normalization
library(sva) # Batch effect removal
library(ggplot2)
library(openxlsx)
library(RColorBrewer)
library(factoextra)
library(ggpubr)
library(clusterProfiler)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "TGCT_tpm+kegg.Rdata")

## Perform correlation analysis
use_data <- as.data.frame(t(TGCT_tpm))
## Create a file to store correlation results
cor_result <- data.frame(matrix(NA, nrow = ncol(use_data), ncol = 3))
colnames(cor_result) <- c("gene", "cor", "pvalue")
rownames(cor_result) <- colnames(use_data)

# Select the target gene
use_gene <- "PRTFDC1"
for (gene_name in colnames(use_data)) {
  
  cor_value <- NA
  p_value <- NA
  
  # Calculate correlation
  cor_value <- cor(as.numeric(use_data[, use_gene]), 
                   as.numeric(use_data[, gene_name]), method = "spearman")
  
  # Calculate p-value for the correlation
  p_value <- cor.test(as.numeric(use_data[, use_gene]), 
                      as.numeric(use_data[, gene_name]), method = "spearman")$p.value
  
  # Save the results
  cor_result[gene_name, 1] <- gene_name
  cor_result[gene_name, 2] <- cor_value
  cor_result[gene_name, 3] <- p_value
}
# Calculate FDR
# The obtained p-values need to be adjusted, using adjusted p-values for more accurate filtering, here using the FDR correction method
cor_result$FDR <- p.adjust(cor_result$pvalue, method = "fdr")
# Filter for significantly correlated genes (adjust criteria as needed)
cor_gene <- cor_result %>% filter(abs(cor) > 0.4, FDR < 0.05)
# Get the list of genes
gene <- cor_gene$gene
## Convert
gene = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

## Remove duplicates
gene <- dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)
gene_df <- data.frame(cor = cor_gene$cor,
                      SYMBOL = cor_gene$gene)
gene_df <- merge(gene_df, gene, by = "SYMBOL")

## geneList trilogy
## 1. Get gene logFC
geneList <- gene_df$cor
## 2. Assign names
names(geneList) = gene_df$SYMBOL
## 3. Sorting is important
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

library(clusterProfiler)

kegg_geneset <- kegg_geneset
### Main GSEA program
y <- GSEA(geneList, TERM2GENE = kegg_geneset)

yd <- as.data.frame(y)

mytheme <- theme(plot.title = element_text(family = "Arial", size = 26, face = "bold", color = "#831A21",
                                           hjust = 0.5,          # Horizontal position of the font
                                           vjust = 0.5,          # Vertical height of the font
                                           angle = 0),           # Font angle
                 plot.subtitle = element_text(family = "Arial", size = 26, face = "bold", color = "black", hjust = 0.5, vjust = 0.5, angle = 0),
                 text = element_text(family = "Arial", size = 26, face = "plain", color = "black"),
                 axis.title = element_text(family = "Arial", size = 16, face = "bold", color = "black"),
                 axis.text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
                 axis.line = element_line(color = "black", size = 1),
                 axis.ticks = element_line(color = "black", size = 1),
                 legend.title = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
                 legend.text = element_text(family = "Arial", size = 14, face = "plain", color = "black"),
                 legend.background = element_blank(),
                 legend.position = 'right') # Custom theme
### View overall distribution
library(ggplot2)
dotplot(y, x = "NES",
        showCategory = 15,
        font.size = 15,
        label_format = 50,
        split = ".sign") + facet_grid(~.sign) + mytheme
head(yd)
msi_hm_symbol <- yd
hm_gse_cut <- msi_hm_symbol %>% filter(pvalue < 0.05 & abs(NES) > 1)
hm_gse_cut_down <- hm_gse_cut[hm_gse_cut$NES < 0,]
hm_gse_cut_up <- hm_gse_cut[hm_gse_cut$NES > 0,]
# You can choose to display the top pathways by NES or other criteria, displaying 10 for each
down_gsea <- hm_gse_cut_down[tail(order(hm_gse_cut_down$NES, decreasing = T), 15),]
up_gsea <- hm_gse_cut_up[head(order(hm_gse_cut_up$NES, decreasing = T), 15),]
diff_gsea <- hm_gse_cut[head(order(abs(hm_gse_cut$NES), decreasing = T), 20),]

p4 <- ggplot(data = up_gsea,
             aes(x = NES, y = reorder(ID, NES))) +  # Reorder pathway based on GeneRatio
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradient(low = '#E42A2A', high = '#14B3FF', name = "p.adjust") +
  scale_size_continuous(range = c(4, 11)) +
  labs(x = "NES",
       y = "",
       title = "Activated",
       size = "Count") +
  theme_bw() +
  mytheme + scale_y_discrete(labels = function(x) str_wrap(x, width = 20))

p4

p4 <- ggplot(data = down_gsea,
             aes(x = NES, y = reorder(ID, -NES))) +  # Reorder pathway based on GeneRatio
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_color_gradient(low = '#E42A2A', high = '#14B3FF', name = "p.adjust") +
  scale_size_continuous(range = c(4, 11)) +
  labs(x = "NES",
       y = "",
       title = "Suppressed",
       size = "Count") +
  theme_bw() +
  mytheme + scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
p4
# Visualization of results
# If you want to display them together
gseaplot2(y,
          geneSetID = c("KEGG_FOCAL_ADHESION",
                        "KEGG_ECM_RECEPTOR_INTERACTION",
                        "KEGG_PATHWAYS_IN_CANCER",
                        "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                        "KEGG_WNT_SIGNALING_PATHWAY"), # Enriched ID number
          title = "Enriched in KEGG", # Title
          color = "red", # GSEA line color
          base_size = 20, # Base font size
          rel_heights = c(2, 0.5, 0.5), # Relative heights of subplots
          subplots = 1:3, # Which subplots to display, e.g., subplots = c(1, 3) # Only the first and third plots
          ES_geom = "line", # Use line or dot for enrichment score
          pvalue_table = T) # Display p-value and other information
gseaplot2(y,
          geneSetID = c("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                        "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
                        "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
                        "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                        "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"), # Enriched ID number
          title = "Enriched in KEGG", # Title
          color = "red", # GSEA line color
          base_size = 20, # Base font size
          rel_heights = c(2, 0.5, 0.5), # Relative heights of subplots
          subplots = 1:3, # Which subplots to display, e.g., subplots = c(1, 3) # Only the first and third plots
          ES_geom = "line", # Use line or dot for enrichment score
          pvalue_table = T) # Display p-value and other information
```

This translation retains all the functionality and comments of your original R script while converting all text into English.