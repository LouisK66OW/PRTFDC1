rm(list = ls())
library(parallel)
library(dplyr)
library(tidyr)
library(PharmacoGx)
library(stringr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)
library(circlize)
library(ComplexHeatmap)
library(limma)
library(org.Hs.eg.db)
Sys.setenv(LANGUAGE = "en") # Display error messages in English
options(stringsAsFactors = FALSE) # Prevent conversion of strings to factors
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## cmap tutorial 1. https://www.jianshu.com/p/1d6623a706a1
##              2. https://www.jianshu.com/p/39a7ed032d8d
##              3. [New CMap database CLUE usage - How to view Detailed List?] https://www.bilibili.com/video/BV1FT411w71G/?share_source=copy_web&vd_source=f24225d34c6a3bcfeeebb989ac24422f

## First, you need to perform differential analysis using signature scores to find the top 150 upregulated and downregulated genes for online cmap analysis
## Load previously analyzed files
load(file = "TCGA-TGCT-mrna_expr+pheno+sur.Rdata")
## Input your gene, then proceed with operations
gene = "PRTFDC1"
## Grouping
## Generate gene grouping data
## Load PCaDB file
ENSEMBL_ID <- bitr(gene, fromType = "SYMBOL", toType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db)[,2]
group <- expr %>% .[ENSEMBL_ID,] %>% t() %>% as.data.frame() 
colnames(group) = "gene"
for (i in 1:length(rownames(group))) {
  if (group$gene[i] >= median(group$gene)) {
    group$group[i] = "High"
  } else {
    group$group[i] = "Low"
  }
}

## Prepare for limma
ann <- group %>% dplyr::select(-gene) 
exprSet <- expr[,rownames(ann)]
group <- factor(ann$group, levels = c("Low", "High")) ## Control group first
## Start limma
### If an error occurs
if (T) {
  # 1. Transpose
  data <- t(exprSet)
  ## 2. Convert NaN values to 0
  data[is.nan(data)] <- 0 
  ## 3. Remove genes with no variance using the var function
  data <- data[, which(apply(data, 2, var) != 0)]
  res.pca <- prcomp(data, scale = TRUE)
}

library(factoextra)
fviz_pca_ind(res.pca, col.ind = group)

### 1. Construct comparison matrix
design <- model.matrix(~group)
### Name the comparison matrix
colnames(design) <- levels(group)
design

### 2. Linear model fitting
fit <- lmFit(exprSet, design)
### 3. Bayesian test
fit2 <- eBayes(fit)
### 4. Output differential analysis results; the number of coef cannot exceed the number of columns in design
### Here, 2 represents the comparison between the second and first columns in design
allDiff = topTable(fit2, adjust = 'fdr', coef = 2, number = Inf) %>% rownames_to_column(var = "ENSEMBLID")
### Add gene names SYMBOL
library(AnnotationDbi)
library(org.Hs.eg.db)
allDiff$symbol <- mapIds(org.Hs.eg.db,
                         keys = allDiff$ENSEMBLID,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
cmapup1 <- allDiff %>% filter(adj.P.Val < 0.05) %>% filter(logFC > 1) ## arrange(desc(logFC)) %>% .[1:150,]
cmapdown1 <- allDiff %>% filter(adj.P.Val < 0.05) %>% filter(logFC < (-1)) ## arrange(logFC) %>% .[1:150,]
cmapupanddown <- rbind(cmapup1, cmapdown1)
cmapup1 <- cmapup1$symbol
cmapdown1 <- cmapdown1$symbol
write.table(cmapup1, file = paste0(gene, "-cmapup.txt"), quote = FALSE, row.names = F, col.names = F)
write.table(cmapdown1, file = paste0(gene, "-cmapdown.txt"), quote = FALSE, row.names = F, col.names = F)
## Use the above 300 genes for calculation using the latest algorithm

## Read in moa; many columns in this gct file are unnecessary, keep only the needed parts
allcmap <- data.table::fread(file = "query_result.gct", skip = 2)[-1,] %>% 
  as.data.frame() %>% 
  dplyr::select(pert_id, pert_iname, pert_type, moa, raw_cs, fdr_q_nlog10, norm_cs, tas) %>% ## Keep only needed parts
  filter(moa != "-666") %>% ## -666 means NA
  filter(pert_type == "trt_cp")
allcmap[,5:8] <- lapply(allcmap[,5:8], as.numeric)
allcmap <- allcmap %>% ## Keep only small molecule drugs
  arrange(norm_cs) %>% ## Sort from low to high
  distinct(pert_id, .keep_all = T) ## Remove duplicate IDs
allcmap <- allcmap[grepl("[A-Z]-", allcmap$pert_iname) == F,]  ## Remove drugs with only codes, optional
## Take the first mechanism if there are multiple
allcmap$moa <- gsub("\\|.*", "", allcmap$moa)
## Calculate FDR
allcmap$FDR <- 10^(-(as.numeric(allcmap$fdr_q_nlog10)))
### Export
write.csv(allcmap, file = "TCGA-TGCT cmap.csv", row.names = F)
## Mimic SCG2 to select the top 20
select_cmap <- allcmap %>% distinct(pert_iname, .keep_all = T) %>% .[1:20,]
## Calculate the count of top 20 inhibitors in 20 cmap
MoAinput <- as.data.frame(table(select_cmap$moa))
colnames(MoAinput) <- c("MoA", "Perturbagen.Count")
for (i in 1:length(rownames(MoAinput))) {
  MoAinput$Perturbagen.Id[i] = paste(select_cmap$pert_id[which(select_cmap$moa == MoAinput$MoA[i])], collapse = ",")
  MoAinput$Name[i] = paste(select_cmap$pert_iname[which(select_cmap$moa == MoAinput$MoA[i])], collapse = ",")
}

# Organize into the input format needed for oncoPrint
PerturbagenID <- unlist(str_split(MoAinput$Name, ","))
names(PerturbagenID) <- unlist(str_split(MoAinput$Perturbagen.Id, ","))
MoAinput <- MoAinput[, c("MoA", "Perturbagen.Id")] %>% split(.$MoA) %>% lapply("[[", 2) %>% 
  lapply(., function(x) unlist(str_split(x, ","))) %>% plyr::ldply(., data.frame)
colnames(MoAinput) <- c("mechanisms of action", "inhibitors")
oncoprintinput <- reshape2::dcast(MoAinput, `mechanisms of action` ~ inhibitors)
rownames(oncoprintinput) <- oncoprintinput$`mechanisms of action`
oncoprintinput <- oncoprintinput[, -1] %>% as.data.frame()
oncoprintinput[!is.na(oncoprintinput)] <- "inhibitor"
oncoprintinput[is.na(oncoprintinput)] <- ""
colnames(oncoprintinput) <- PerturbagenID[colnames(oncoprintinput)]
oncoprintinput <- oncoprintinput[, order(colnames(oncoprintinput))]
select_cmap <- select_cmap %>% arrange(pert_iname)
oncoprintinput <- oncoprintinput[, select_cmap$pert_iname]
select_cmap <- select_cmap %>% arrange(norm_cs)
oncoprintinput <- oncoprintinput[, select_cmap$pert_iname]

### Plot

alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "white", col = "grey")),
  # dots
  inhibitor = function(x, y, w, h) 
    grid.points(x, y, pch = 16, size = unit(1.2, "char"))
)
ha_coldata <- colSums(apply(oncoprintinput, 2, function(x) x == "inhibitor") + 0) %>% as.numeric()
ha_rowdata <- rowSums(apply(oncoprintinput, 2, function(x) x == "inhibitor") + 0) %>% as.numeric() 
top_ha <- HeatmapAnnotation(inhibitors = anno_barplot(ha_coldata, axis = F, border = F, 
                                                      gp = gpar(fill = "grey"),
                                                      bar_width = 1),
                            annotation_name_side = "left",
                            annotation_name_rot = 90)
right_ha <- rowAnnotation(Count = anno_barplot(ha_rowdata, axis = F, border = F, 
                                               gp = gpar(fill = "grey"),
                                               bar_width = 1, width = unit(1.5, "cm")),
                          annotation_name_side = "top",
                          annotation_name_rot = 0)

pdf("MoA.pdf", width = 7, height = 5, onefile = FALSE)
#png("MoA.png", width = 480, height = 480)
oncoPrint(oncoprintinput, alter_fun = alter_fun, 
          show_column_names = TRUE, column_names_side = "top",
          column_order = 1:ncol(oncoprintinput), 
          top_annotation = top_ha,
          right_annotation = right_ha,
          show_pct = FALSE, show_heatmap_legend = F, alter_fun_is_vectorized = FALSE)
decorate_annotation("inhibitors", {
  grid.text("Mechanism of Action", unit(1, "npc") + unit(3, "mm"), just = "left")})
dev.off()

# Plot bubble chart
library(gcookbook)
select_cmap$`-Normalized Connectivity Scores` <- (-(select_cmap$norm_cs))
pp = ggplot(select_cmap, aes(x = `-Normalized Connectivity Scores`, y = reorder(pert_iname, -norm_cs)))

pp + ggtitle("TCGA-TGCT") +
  geom_point(aes(color = tas, size = `-Normalized Connectivity Scores`)) +
  geom_segment(aes(x = signif(`-Normalized Connectivity Scores`, 2)[20], xend = (`-Normalized Connectivity Scores` - 0.006), y = pert_iname, yend = pert_iname), size = 1, color = "black") + 
  theme_bw() +
  ylab(NULL) +
  xlab("-NCS") +
  scale_colour_gradient(low = "#D77071", high = "#222f3e") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20), position = "left") +
  labs(colour = "tas", size = "-NCS") +
  theme(axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(color = "black", size = 14, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 15, face = 2), legend.position = "right"
  ) #+scale_x_reverse()
ggsave("TCGA-TGCT bubble.pdf", width = 7, height = 7)
ggsave("TCGA-TGCT bubble.png", width = 7, height = 7)
dev.off()
