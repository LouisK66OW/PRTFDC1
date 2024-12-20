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

load(file = "PCA_batch.Rdata")

############## PCA

# rt<-PCA3218_10783_96_batch_After
# rt<-PCA3218_10783_96_batch_Before
# rt<-PCA3218_10783_97_batch_After
# rt<-PCA3218_10783_97_batch_Before
rt <- rt %>% 
  column_to_rownames("sample")
data = rt[, c(2:ncol(rt))]
Type = rt[, 1]  # Extract grouping information
var = colnames(rt)[1]
# PCA analysis
data.pca = prcomp(data, scale. = TRUE)
library(ggplot2)
library(factoextra)
library(FactoMineR)
fviz_pca_ind(data.pca,
             geom.ind = "point", # Show only points
             pointsize = 5, # Size of points
             habillage = rt$group, # Shape of points
             fill.ind = rt$group, # Group color
             palette = c("#E93639", "#54D9ED"),
             addEllipses = TRUE, ellipse.level = 0.95,  # Add confidence ellipses
             legend.title = "Groups", # Legend title
             title = "PCA-GSE3218-GSE10783",
             subtitle = "GPL97") +
  theme_bw() + # Beautify with ggplot2
  theme(plot.title = element_text(family = "Arial", size = 20, face = "bold", color = "black",
                                  hjust = 0.5,          # Horizontal position of font
                                  vjust = 0.5,          # Vertical height of font
                                  angle = 0),  # Angle of font tilt
        plot.subtitle = element_text(family = "Arial", size = 20, face = "bold", color = "black", hjust = 0.5, vjust = 0.5, angle = 0),  # Angle of font tilt
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
