Here is the translated version of your R code with comments and text in English:
  
  ```r
rm(list = ls())
## Load R packages
options(stringsAsFactors = F)
library(data.table)
library(stringr)
library(Biobase)
library(BiocGenerics)
library(GEOquery)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(devtools) # Install R packages from GitHub
library(limma) # Normalization
library(sva) # Batch correction
library(ggplot2)
library(openxlsx)
library(scatterplot3d)
library(pheatmap)
library(ggplot2)
library(plot3D)
library(RColorBrewer)
library(ggplot2)
library(factoextra)
library(FactoMineR)
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load(file = "expr+diff+alldiff.Rdata")
############## PCA Analysis
# rt <- PCA3218_10783_96group
# rt <- PCA3218_10783_97group
rt <- rt %>% 
  column_to_rownames("sample")
data = rt[, 2:ncol(rt)]
Type = rt[, 1]  # Extract group information
var = colnames(rt)[1]
# PCA analysis
data.pca = prcomp(data, scale. = TRUE)
library(ggplot2)
library(factoextra)
library(FactoMineR)
fviz_pca_ind(data.pca,
             geom.ind = "point", # show points only
             pointsize = 5, # size of points
             habillage = rt$group, # shape of points
             fill.ind = rt$group, # colors by group
             palette = c("#E93639", "#54D9ED"),
             addEllipses = TRUE, ellipse.level = 0.95, # add confidence ellipses
             legend.title = "Groups",
             title = "PCA-GSE3218-GSE10783",
             subtitle = "GPL97") +
  theme_bw() + # Beautify with ggplot2
  theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", color = "black",
                                  hjust = 0.5, vjust = 0.5, angle = 0),
        plot.subtitle = element_text(family = "Arial", size = 20, face = "bold", color = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        axis.title = element_text(family = "Arial", size = 18, face = "bold", color = "black"),
        axis.text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 1),
        legend.title = element_text(family = "Arial", size = 18, face = "plain", color = "black"),
        legend.text = element_text(family = "Arial", size = 18, face = "plain", color = "black"),
        legend.background = element_blank(),
        legend.position = "right"
  )

########### Differential Analysis
group <- c(rep("cancer", 101), rep("normal", 6), rep("cancer", 34))
group <- as.character(group)
### Convert group to vector and specify levels order
### In levels, put the control group first
group <- factor(group, levels = c("normal", "cancer"))
table(group)

### Construct comparison matrix
design <- model.matrix(~group)
### Name the comparison matrix
colnames(design) <- levels(group)
design

### 2. Linear model fitting
fit <- lmFit(expr_merge, design)
### 3. Bayesian testing
fit2 <- eBayes(fit)
### 4. Output differential analysis results, where coef number should not exceed the number of columns in design
### Here, 2 represents the comparison between the second and first columns in design
allDiff = topTable(fit2, adjust = 'fdr', coef = 2, number = Inf)

###################################################################################
### Define differentially expressed genes: fold change (logFC) > 1.5, adjusted p-value < 0.05

colnames(allDiff)

library(dplyr)
diffgene <- allDiff %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC) > 1.5)

############## Volcano Plot
# Get gene list
allDiff <- allDiff_96 %>% 
  ## Convert row names to column names, add a new column
  rownames_to_column("ENSEMBL")
## Conversion
#### Use ensemble to match symbol
library(org.Hs.eg.db)
k = keys(org.Hs.eg.db, keytype = "SYMBOL")
head(k, 5)
list = AnnotationDbi::select(org.Hs.eg.db, keys = k, columns = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")
allDiff_1 <- merge(list, allDiff, by = "ENSEMBL")
allDiff_1 <- allDiff_1 %>% 
  ## Remove redundant information
  dplyr::select(-ENSEMBL) %>% 
  ## Sort expression levels in descending order
  arrange(desc(logFC)) %>% 
  ## Remove duplicates, keep the first symbol
  distinct(SYMBOL, .keep_all = T) %>% 
  ## Convert column names to row names
  column_to_rownames("SYMBOL")
allDiff_1 <- allDiff_1 %>% 
  arrange(logFC)
data <- allDiff_1
ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  # Horizontal reference line:
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  # Vertical reference lines:
  geom_vline(xintercept = c(-1.2, 1.2), linetype = "dashed", color = "#999999") +
  # Scatter plot:
  geom_point(aes(size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  # Specify color gradient mode:
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  # Specify size gradient mode:
  scale_size_continuous(range = c(1, 3)) +
  # Theme adjustment:
  theme_bw() +
  # Adjust theme and legend position:
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.7),
        legend.justification = c(0, 1)
  ) +
  # Set some legends not to display:
  guides(col = guide_colourbar(title = "-Log10(adj.P-value)"),
         size = "none") +
  # Add labels:
  geom_text(aes(label = label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust = 1) +
  # Modify axes:
  xlab("Log2(Fold Change)") +
  ylab("-Log10(adj.P-value)") +
  
  ggtitle("GSE3218-GSE10783",
          subtitle = "GPL96") +
  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", color = "black",
                                  hjust = 0.5, vjust = 0.5, angle = 0),
        plot.subtitle = element_text(family = "Arial", size = 20, face = "bold", color = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        axis.title = element_text(family = "Arial", size = 18, face = "bold", color = "black"),
        axis.text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 1),
        legend.title = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        legend.text = element_text(family = "Arial", size = 16, face = "plain", color = "black"),
        legend.background = element_blank(),
        legend.position = "right"
  )

############## Functional Enrichment Analysis

# Get gene list
diffgene_96 <- diffgene_96 %>% 
  ## Convert row names to column names, add a new column
  rownames_to_column("ENSEMBL")
diffgene_97 <- diffgene_97 %>% 
  ## Convert row names to column names, add a new column
  rownames_to_column("ENSEMBL")
diffgene <- rbind(diffgene_96, diffgene_97)
res1 <- diffgene
## Filter differentially expressed genes
# First, sort the table by padj value in ascending order, continue sorting by log2FC in descending order under the same padj value
res1 <- res1[order(res1$P.Value, res1$logFC, decreasing = c(FALSE, TRUE)), ]
# Mark others as none, representing non-differentially expressed genes
res1[which(res1$logFC >= 1.5 & res1$P.Value < 0.05), 'sig'] <- 'up'
res1[which(res1$logFC <= -1.5 & res1$P.Value < 0.05), 'sig'] <- 'down'
res1[which(abs(res1$logFC) <= 1.5 | res1$P.Value >= 0.05), 'sig'] <- 'none'
# Output the selected differential gene master table
res1_select <- subset(res1, sig %in% c('up', 'down'))
# Output up and down separately
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

diffgene <- res1_up
# diffgene <- res1_down

## Conversion
#### Use ensemble to match symbol
library(org.Hs.eg.db)
k = keys(org.Hs.eg.db, keytype = "SYMBOL")
head(k, 5)
list = AnnotationDbi::select(org.Hs.eg.db, keys = k, columns = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")
diffgene <- merge(list, diffgene, by = "ENSEMBL")

diffgene <- diffgene %>% 
  ## Remove redundant information
  dplyr::select(-ENSEMBL) %>% 
  ## Sort expression levels in descending order
  arrange(desc(logFC)) %>% 
  ## Remove duplicates, keep the first symbol
  distinct(SYMBOL, .keep_all = T) %>% 
  ## Convert column names to row names
  column_to_rownames("SYMBOL")
diff <- diffgene
gene <- rownames(diff)
# Convert the entire differential gene set as an example:
## Conversion of up-regulated or down-regulated genes is the same, not repeated here
diff_entrez <- bitr(gene,
                    fromType = "SYMBOL", # Current ID type
                    toType = "ENTREZID", # ID type to convert
                    OrgDb = "org.Hs.eg.db")
# KEGG enrichment analysis (hypergeometric distribution test):
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa", # Species Homo sapiens
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500)

# Convert ENTREZ back to symbol:
KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
View(KEGG_diff@result)
colnames(KEGG_diff@result) # Explain the meanings of different column names

# Calculate Rich Factor:
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# Calculate Fold Enrichment:
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]
# Extract the KEGG enrichment analysis result table:
KEGG_result <- KEGG_diff2@result
# Example with Top 20 results from the enrichment table:
KEGG_top20 <- KEGG_result[1:20,]
## Calculate GeneRatio, first split GeneRatio into two fields, calculate the ratio
KEGG_top20  <- separate(KEGG_top20, col="GeneRatio", into=c("count", "pathway_count"), sep="[/]", remove=T)
KEGG_top20$count <- as.numeric(KEGG_top20$count)
KEGG_top20$pathway_count <- as.numeric(KEGG_top20$pathway_count)
KEGG_top20$GeneRatio <- round(KEGG_top20$count / KEGG_top20$pathway_count, 2)

# Specify plotting order (convert to factor):
KEGG_top20$pathway <- factor(KEGG_top20$Description, levels = rev(KEGG_top20$Description))

mytheme <- theme(plot.title=element_text(family = "Arial", size=22, face="bold", color="#831A21",
                                         hjust = 0.5,          # Horizontal position of the font
                                         vjust = 0.5,          # Vertical height of the font
                                         angle = 0),           # Angle of the font
                 plot.subtitle = element_text(family = "Arial", size=18, face="bold", color="black", hjust = 0.5, vjust = 0.5, angle = 0),
                 text=element_text(family = "Arial", size=14, face="plain", color="black"),
                 axis.title=element_text(family = "Arial", size=16, face="bold", color="black"),
                 axis.text = element_text(family = "Arial", size=16, face="plain", color="black"),
                 axis.line = element_line(color = "black", size = 1),
                 axis.ticks = element_line(color = "black", size = 1),
                 legend.title = element_text(family = "Arial", size=16, face="plain", color="black"),
                 legend.text = element_text(family = "Arial", size=14, face="plain", color="black"),
                 legend.background = element_blank(),
                 legend.position='right') # Custom theme

p4 <- ggplot(data = KEGG_top20,
             aes(x = GeneRatio, y = reorder(Description, Count))) +  # Reorder pathway based on GeneRatio
  geom_point(aes(size = Count, color = -log10(pvalue))) +
  scale_colour_gradient(low="yellow", high="red") +
  scale_size_continuous(range = c(4, 11)) +
  labs(x = "GeneRatio",
       y = "",
       title = "up",
       size = "Count") +
  theme_bw() +
  mytheme
p4

p4 + scale_y_discrete(labels=function(x) str_wrap(x, width=40))

