R
library(openxlsx)
library(maftools)
library(xlsx)
library(tidyr)
library(dplyr)
library(limma)
library(readr)
library(tibble)
library(UCSCXenaTools)
library(survminer) # Load package
library(survival) # Load package
library(maftools)
library(TCGAbiolinks)
library(ComplexHeatmap)
library(egg)
library(cowplot)
library(RColorBrewer)
library(aplot)
library(clusterProfiler)
library(org.Hs.eg.db)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##Enter your cancer type
cas <- "TGCT"
##Enter your gene
gene <- "PRTFDC1"
try <-  function(ca){
  library(openxlsx)
  library(maftools)
  library(xlsx)
  library(tidyr)
  library(dplyr)
  library(limma)
  library(readr)
  library(tibble)
  library(UCSCXenaTools)
  library(survminer) # Load package
  library(survival) # Load package
  library(maftools)
  library(TCGAbiolinks)
  library(ComplexHeatmap)
  library(egg)
  library(cowplot)
  library(RColorBrewer)
  library(aplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  ##Enter your gene
  genedds <- gene
  gene <- bitr(gene, fromType ="SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db, drop = F)[1,2]
  ##Enter your cancer type
  # ca <- ""
  gdcca <- unique(XenaData$XenaCohorts[grepl("GDC", XenaData$XenaCohorts) & grepl(ca, XenaData$XenaCohorts)])
  #XenaData <- XenaData %>% filter(XenaCohorts==gdcca)
  ###Check what's available
  #view(XenaData)
  ##If not available, download
  ## Download FPKM expression matrix data
  #eset <- XenaGenerate(subset = XenaCohorts == gdcca) %>% 
  #  XenaFilter(filterDatasets = paste0(paste0("TCGA-", ca), ".htseq_fpkm.tsv")) %>% 
  #  XenaQuery() %>%
  #  XenaDownload(destdir = getwd()) %>% 
  #  XenaPrepare()
  ##Download clinical data
  #sur <- XenaGenerate(subset = XenaCohorts == gdcca) %>% 
  #  XenaFilter(filterDatasets = paste0(paste0("TCGA-", ca), ".survival.tsv")) %>% 
  #  XenaQuery() %>%
  #  XenaDownload(destdir = getwd()) %>% 
  #  XenaPrepare()
  #
  ##Download phenotype data
  #pheno <- XenaGenerate(subset = XenaCohorts == gdcca) %>% 
  #  XenaFilter(filterDatasets = paste0(paste0("TCGA-", ca), ".GDC_phenotype.tsv")) %>% 
  #  XenaQuery() %>%
  #  XenaDownload(destdir = getwd()) %>% 
  #  XenaPrepare()
  ##Download CNV
  #CNV <- XenaGenerate(subset = XenaCohorts == gdcca) %>% 
  #  XenaFilter(filterDatasets = paste0(paste0("TCGA-", ca), ".masked_cnv.tsv")) %>% 
  #  XenaQuery() %>%
  #  XenaDownload(destdir = getwd()) %>% 
  #  XenaPrepare()
  #
  ###Start downloading multi-omics data, sometimes a VPN is needed
  ##query <- GDCquery(
  ##  project = paste0("TCGA-", ca), 
  ##  data.category = "Simple Nucleotide Variation",
  ##  data.type = "Masked Somatic Mutation",
  ##  cases = "open"
  ##)
  ##GDCdownload(query, method = "api", files.per.chunk = NULL)
  ##GDCprepare(query, save = T, save.filename = paste0(paste0("TCGA-", ca), "_SNP.Rdata"))
  
  ##First unzip into loose files
  #setwd(dir = ca)
  #files <- list.files(pattern = '*.gz', recursive = TRUE)
  #data <- data.frame()
  #for (file in files) {
  #  mut <- read.delim(file, skip = 7, header = T, fill = TRUE, sep = "\t")
  #  data <- rbind(data, mut)
  #}
  #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  ###If already downloaded, load the data
  #eset <- data.table::fread(file = paste0(paste0("TCGA-", ca), ".htseq_fpkm.tsv.gz"))
  #sur <- data.table::fread(paste0(paste0("TCGA-", ca), ".survival.tsv"))
  #pheno <- data.table::fread(paste0(paste0("TCGA-", ca), ".GDC_phenotype.tsv.gz"))
  ###Start
  #test <- pheno %>% filter(sample_type.samples == "Primary Tumor") %>% distinct(patient_id, .keep_all = T)
  #table(test$submitter_id.samples %in% colnames(eset))
  #expr <- eset %>% column_to_rownames(var = colnames(eset)[1]) %>% .[, colnames(.) %in% test$submitter_id.samples]
  #rownames(expr) <- substring(rownames(expr), 1, 15)
  #
  ###Save
  #save(expr, sur, data, pheno, file = paste0(paste0("TCGA-", ca), " SNP.Rdata"))
  load(file = paste0(paste0("TCGA-", ca), " SNP.Rdata"))
  
  
  ##Colors: high #EE2C2C Low #87CEEB
  "#EE2C2C"
  "#87CEEB"
  
  ##Try survival curve
  #names(sur)[which(colnames(sur) == "sample")] <- "sample_id"
  #sur <- inner_join(sur, group, by = "sample_id")
  #fit <- surv_fit(Surv(OS.time, OS) ~ group,
  #                data = sur)
  #ggsurvplot(fit, data = sur,
  #           conf.int = TRUE,  # Add confidence interval
  #           risk.table = TRUE, # Draw cumulative risk curve
  #           pval = TRUE, # Add P-value
  #           palette = c("#EE2C2C", "#87CEEB")) 
  
  
  
  
  
  ##Start from FMGS
  
  
  
  ##aow is to ensure the control group is on the left, low is automatically on the right, modify in AI after the graph is generated
  
  maf <- data[data$Variant_Classification != "Silent", ]   # Filter silent mutation
  mut <- pmin(table(maf$Hugo_Symbol, maf$Tumor_Sample_Barcode), 1)
  mut <- as.data.frame(mut)
  mut <- data.table::dcast(mut, formula = Var1 ~ Var2)
  rownames(mut) <- mut$Var1; mut <- dplyr::select(mut, -Var1)
  
  
  
  # Mutation main region grouping annotation
  ##Start calculating gene and its group
  group <- expr[gene, colnames(expr) %in% substring(colnames(mut), 1, 16)[!duplicated(substring(colnames(mut), 1, 16))]] %>% 
    t() %>%
    as.data.frame() %>% 
    rownames_to_column(var = "sample_id") 
  colnames(group)[2] <- "gene"
  for (i in 1:length(rownames(group))) {if(group$gene[i] >= median(group$gene)){group$group[i] = "High"}
    else{group$group[i] = "Low"}
    
  }
  
  mut_gene <- mut[, !duplicated(substring(colnames(mut), 1, 16))]
  colnames(mut_gene ) <- substring(colnames(mut_gene), 1, 16)  
  mut_gene <- mut_gene[, group$sample_id]
  dataa <- data[substring(data$Tumor_Sample_Barcode, 1, 16) %in% group$sample_id,] %>% .[.$Variant_Classification != "Silent", ] 
  
  
  ##Top 10 genes
  clin_for_maf <- dataa[, "Tumor_Sample_Barcode", drop = F] %>% distinct(Tumor_Sample_Barcode, .keep_all = T) %>% dplyr::mutate(`submitter_id.samples` = substring(Tumor_Sample_Barcode, 1, 16))
  for (i in 1:nrow(clin_for_maf)) {if(substring(clin_for_maf$Tumor_Sample_Barcode[i], 1, 16) %in% group$sample_id[which(group$group == "High")] == T){
    clin_for_maf$Group[i] = "High"}
    else{clin_for_maf$Group[i] = "aow"}}
  pheno <- pheno[pheno$submitter_id.samples %in% clin_for_maf$submitter_id.samples,]
  #Rename stage
  pheno$new_stage = ifelse(pheno$clinical_stage %in% c('Stage I', 'Stage IA', 'Stage IB', "Stage IS"), 's1',
                           ifelse(pheno$clinical_stage %in% c('Stage II', 'Stage IIA', 'Stage IIB', 'Stage IIC'), 's2',
                                  ifelse(pheno$clinical_stage %in% c('Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC'), 's3',
                                         ifelse(pheno$clinical_stage %in% c('Stage IV', 'Stage IVA', 'Stage IVB', 'Stage IVC'), 's4', 'NA'
                                         )  ) ) )
  clin_for_maf <- inner_join(clin_for_maf, pheno[, c("submitter_id.samples", "new_stage", "clinical_T", "clinical_N", "clinical_M")], by = "submitter_id.samples") %>% dplyr::rename(Stage = new_stage, T = clinical_T, N = clinical_N, M = clinical_M)
  ###Colors
  vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
  names(vc_cols) = c(
    'Frame_Shift_Del',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Multi_Hit',
    'Frame_Shift_Ins',
    'In_Frame_Ins',
    'Splice_Site',
    'In_Frame_Del'
  )
  #Color and variant types must match in number
  Scolors = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
  names(Scolors) = unique(clin_for_maf$Stage)
  
  Tc <- RColorBrewer::brewer.pal(n = length(unique(clin_for_maf$T)), name = 'Reds')
  names(Tc) = unique(clin_for_maf$T)
  
  Nc <- RColorBrewer::brewer.pal(n = length(unique(clin_for_maf$N)), name = 'Blues')
  names(Nc) = unique(clin_for_maf$N)
  
  Mc <- RColorBrewer::brewer.pal(n = length(unique(clin_for_maf$M)), name = 'Greens')
  names(Mc) = unique(clin_for_maf$M)
  
  colors = list(Group = c(aow = '#87CEEB', High = '#EE2C2C'), Stage = Scolors, T = Tc, N = Nc, M = Mc)
  ##
  maf <- read.maf(dataa, clinicalData = clin_for_maf)
  genes <- maf@gene.summary$Hugo_Symbol[1:10]
  
  maftools::oncoplot(maf, genes = genes, colors = vc_cols, clinicalFeatures = c("Group", "Stage", "T", "N", "M"), sortByAnnotation = TRUE, showTitle = F, groupAnnotationBySize = F, drawColBar = F, anno_height = 2, annotationColor = colors, writeMatrix = T)#list(Group = c(aow = '#87CEEB', High = '#EE2C2C')))
  export::graph2pdf(file = paste0(genedds, paste0(paste0("TCGA-", ca), " maftools.pdf")), width = 9, height = 6)
  
  High_maf <- dataa[dataa$Tumor_Sample_Barcode %in% clin_for_maf$Tumor_Sample_Barcode[which(clin_for_maf$Group == "High")],]
  Low_maf <- dataa[dataa$Tumor_Sample_Barcode %in% clin_for_maf$Tumor_Sample_Barcode[which(clin_for_maf$Group == "aow")],]
  High_maf <- read.maf(High_maf)
  Low_maf <- read.maf(Low_maf)
  comp <- mafCompare(m1 = High_maf,
                     m2 = Low_maf,
                     m1Name = "High group",
                     m2Name = "Low group",
                     minMut = 0)
  ca_result <- as.data.frame(comp$results)
  mut_gene <- mut_gene[genes,]
  group$group <- factor(group$group, levels = c("Low", "High"), ordered = TRUE)
  group$group <- gsub("Low", "aow", group$group)
  ##Heatmap
  result = ca_result[ca_result$Hugo_Symbol %in% genes, 1:4]
  result$High = round(result$High / length(clin_for_maf$Group == "High"), 2)
  result$Low = round(result$Low / length(clin_for_maf$Group == "aow"), 2)
  dev.off()
  A <- result[, c("Hugo_Symbol", "High", "Low")] %>% pivot_longer(cols = !Hugo_Symbol,
                                                                  names_to = "Cohort",
                                                                  values_to = "frequency")
  mutgene <- character()
  
  
  mutgene <- character()
  for (i in 1:length(genes)) { 
    g = genes[i]
    if (ca_result$pval[ca_result == g] < 0.001) {
      mutgene[i] = paste0("*** ", g)
    } else if (ca_result$pval[ca_result == g] < 0.01) {
      mutgene[i] = paste0("** ", g)
    } else if (ca_result$pval[ca_result == g] < 0.05) {
      mutgene[i] = paste0("* ", g)
    } else if (ca_result$pval[ca_result == g] > 0.05) {
      mutgene[i] = g
    }
  }
  
  for (i in 1:length(A$Hugo_Symbol)) {
    g = A$Hugo_Symbol[i]
    if (ca_result$pval[ca_result == g] < 0.001) {
      A$Hugo_Symbol[i] = paste0("*** ", g)
    } else if (ca_result$pval[ca_result == g] < 0.01) {
      A$Hugo_Symbol[i] = paste0("** ", g)
    } else if (ca_result$pval[ca_result == g] < 0.05) {
      A$Hugo_Symbol[i] = paste0("* ", g)
    } else if (ca_result$pval[ca_result == g] > 0.05) {
      A$Hugo_Symbol[i] = g
    }
  }
  
  A$Hugo_Symbol <- factor(A$Hugo_Symbol, levels = rev(mutgene))
  
  A$Cohort <- factor(A$Cohort, levels = c("Low", "High"))
  p1 <- ggplot(A, aes(x = Cohort, y = Hugo_Symbol, fill = `frequency`)) +
    geom_raster() +
    geom_tile(lwd = 1, linetype = 1, colour = "black") +
    theme_minimal() +
    scale_fill_gradient(limits = c(0, 1), low = 'white', high = 'red') +
    geom_text(aes(label = sprintf("%.3f", `frequency`)), col = 'black', cex = 4) +
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title = element_blank()) + 
    scale_y_discrete(position = "left")
  
  p1
  
  group_heatmap <- ggplot(data = A) + 
    geom_tile(aes(x = Cohort, y = 1, fill = Cohort), lwd = 1, linetype = 1, colour = "black") +
    scale_fill_manual(values = c(Low = '#87CEEB', High = '#EE2C2C')) +
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), axis.line = element_blank(), 
          axis.ticks = element_blank(), axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.position = "right") 
  
  group_heatmap
  
  p2 <- p1 %>% insert_top(group_heatmap, height = .05)
  print(p2)
  export::graph2pdf(file = paste0(genedds, paste0(paste0("TCGA-", ca), "heatmap_FMGs.pdf")), width = 3.2, height = 5)
  dev.off()
  
  #### Mutation boxplot
  gene_matrix <- mut_gene
  gene_matrix[gene_matrix == ""] = "Wild"
  gene_matrix[gene_matrix == 0] = "Wild"
  gene_matrix[gene_matrix != "Wild"] = "Mut"
  gene_matrix = data.frame(t(gene_matrix))
  gene_matrix = rownames_to_column(gene_matrix, var = "sample_id")
  gene_matrix = inner_join(group[, c("gene", "sample_id")], gene_matrix, by = "sample_id")
  gene_matrix = column_to_rownames(gene_matrix, var = "sample_id")
  gene_matrix = gather(gene_matrix, "mutgene", "Group", 2:ncol(gene_matrix))
  colnames(gene_matrix)[1] = paste0(genedds, " Expression")
  gene_matrix$Group = factor(gene_matrix$Group, c("Wild", "Mut"))
  gene_matrix$mutgene = factor(gene_matrix$mutgene, genes) ## Fix gene order
  p <- ggboxplot(gene_matrix, x = "Group", y = paste0(genedds, " Expression"), color = "Group", 
                 palette = "lancet", # Use "jco" for plot, or change to "lancet"
                 add = "jitter", facet.by = "mutgene", short.panel.labs = T,
                 size = 0.3, width = 0.5) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 3) # R defaults to method = "wilcox.test" # Add p-value, running this line will get the following plot
  
  print(p)
  export::graph2pdf(file = paste0(genedds, paste0(paste0(" TCGA-", ca), "batch mut genes.pdf")), width = 6, height = 6)
  dev.off()
  
  #### CNV
  ### Filter mRNA
  # load(file="gtf_df42.Rdata")
  # mrna_gtf_df42 <- gtf_df %>% filter(gene_type == "protein_coding") %>% filter(type == "gene") %>% distinct(gene_id, .keep_all = T) %>% dplyr::select(gene_id, gene_type, gene_name)
  # mrna_gtf_df42$gene_id <- substring(mrna_gtf_df42$gene_id, 1, 15)
  # TGCT_mran_expr <- expr %>% rownames_to_column(var = "gene_id") %>% merge(mrna_gtf_df42[,-2], ., by = "gene_id")
  
  ## Organize CNV files
  # TGCT_CNV <- data.table::fread(file="TCGA-TGCT.masked_cnv.tsv.gz")
  # TGCT_High_CNV <- TGCT_CNV[TGCT_CNV$sample %in% group$sample_id[group$group == "High"]]
  # TGCT_Low_CNV <- TGCT_CNV[TGCT_CNV$sample %in% group$sample_id[group$group == "Low"]]
  # TGCT_All_CNV <- TGCT_CNV[TGCT_CNV$sample %in% group$sample_id]
  # write.table(TGCT_High_CNV, file="TGCT_High_CNV.seg.txt", sep="\t", row.names=F, quote = F)
  # write.table(TGCT_Low_CNV, file="TGCT_Low_CNV.txt", sep="\t", row.names=F, quote = F)
  # write.table(TGCT_All_CNV, file="TGCT_All_CNV.txt", sep="\t", row.names=F, quote = F)
  #### Proceed to GISTIC analysis
  
  
  ## Load GISTIC results
  TGCT.gistic = readGistic(gisticAllLesionsFile = "564158/all_lesions.conf_90.txt", ## It might be 90, check yourself
                           gisticAmpGenesFile = "564158/amp_genes.conf_90.txt", 
                           gisticDelGenesFile = "564158/del_genes.conf_90.txt", 
                           gisticScoresFile = "564158/scores.gistic")
  library(maftools)
  ## 'aow' is to ensure the control group is on the left, 'low' will automatically be on the right, adjust with AI after plotting
  ## Continue to import TMB values
  ## Copy number GISTIC2.0 results (obtained using SNP segment file on GenePattern)
  cna.region <- read.table("564158/all_lesions.conf_90.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
  cna.gene <- read.table("564158/all_thresholded.by_genes.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
  ## Group annotation
  my_ann <- group %>% dplyr::arrange(gene) %>% dplyr::select(sample_id, group) %>% column_to_rownames(var = "sample_id")
  clust.col <- c("#87CEEB", "#EE2C2C")
  manual_order = c("aow", "High")
  my_annotation = HeatmapAnnotation(df = my_ann, 
                                    col = list(group = c("aow" = clust.col[1], "High" = clust.col[2])))
  
  ### Up to arm level
  cna <- cna.region[1:(nrow(cna.region)/2), c(1, 8, 9:(ncol(cna.region)-1))] %>% distinct(Descriptor, .keep_all = T) # Select valid columns
  rownames(cna) <- paste0(gsub(" ", "", cna$Descriptor), "-", substr(rownames(cna), 1, 3)) # Rename rows to identify amplification and deletion sites
  cna.modified <- cna[1:nrow(cna), 3:ncol(cna)] %>% .[, colnames(mut_gene)]
  test_cna <- cna.modified
  test_cna[test_cna != 0] <- "CNV"
  test_cna[test_cna == 0] <- "Not CNV"
  onco.input2 <- cna.modified
  tmp1 <- onco.input2[grepl("Amp", rownames(cna.modified)),] 
  tmp1[tmp1 == 1] <- "Gain" # Values greater than 0 are Gain
  tmp1[tmp1 == 2] <- "Gain"
  tmp1[tmp1 == 0] <- ""
  
  tmp2 <- onco.input2[grepl("Del", rownames(cna.modified)),] # Last 9 are deletions
  tmp2[tmp2 == 1] <- "Loss"
  tmp2[tmp2 == 2] <- "Loss"
  tmp2[tmp2 == 0] <- ""
  onco.input2 <- rbind.data.frame(tmp1, tmp2)
  target <- group %>% dplyr::select(-group) %>% column_to_rownames(var = "sample_id") %>% .[colnames(mut_gene), , drop = FALSE] %>% t() %>% as.data.frame() 
  test2 <- rbind(target, test_cna) %>% t() %>% as.data.frame()
  color = c("#B2DFEE", "#FFFF00")
  test2$gene <- as.numeric(test2$gene)
  for (i in 1:length(rownames(test_cna))) { 
    p = i + 1
    
    # Check if test2[, p] is all "CNV"
    if(all(test2[, p] == "CNV")) {
      test_cna$`p value`[i] = 1
    } else {
      # Convert test2[, p] to a factor with only two levels
      group_factor <- factor(test2[, p], levels = unique(test2[, p]))
      
      # Perform Wilcoxon rank-sum test
      test_cna$`p value`[i] = wilcox.test(test2$gene ~ group_factor, data = test2, var.equal = TRUE)[[3]]
    }
    
    # Count the number of "Gain" and "Loss"
    test_cna$Gain[i] = sum(onco.input2[i, ] == "Gain", na.rm = TRUE)
    test_cna$Loss[i] = sum(onco.input2[i, ] == "Loss", na.rm = TRUE)
    
    # Count the number of "CNV"
    test_cna$CNV[i] = table(test_cna[i, ] == "CNV")[2]
    
    print(i)
  }
  view(test_cna[, (ncol(test_cna)-3):ncol(test_cna)])
  
  ## Select top ten for each
  table <- test_cna[which(test_cna$`p value` != 1), (ncol(test_cna)-3):ncol(test_cna)] %>% dplyr::arrange(`p value`) %>% dplyr::arrange(desc(Gain)) 
  gain10 <- rownames(table)[1:10]
  table <- table %>% dplyr::arrange(`p value`) %>% dplyr::arrange(desc(Loss)) 
  loss10 <- rownames(table)[1:10]
  TGCT.gistic_arm <- c(gain10, loss10) ##
  ## Prepare for plotting
  blue <- "#B2DFEE" #"#5bc0eb"
  red <- "orangered" #"#FFFF00"
  alter_fun2 = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#dcddde", col = "#dcddde"))
    },
    Gain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = red, col = red)) 
    },
    Loss = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = blue, col = blue)) 
    }
  )
  col2 = c("Gain" = red,
           "Loss" = blue)
  # Determine the order of display
  lesion.order <- TGCT.gistic_arm
  tmp <- as.data.frame(t(cna.modified[lesion.order, rownames(my_ann), ]))
  tmp[tmp > 0] <- 1 # Change all values greater than 1 to 1 for calculating mutation frequency
  tmp$group <- as.character(my_ann[rownames(tmp), "group"])
  pct <- NULL
  for (i in lesion.order) {
    tmp1 <- tmp[, c(i, "group")]
    tmp1 <- as.data.frame.array(table(tmp1[, 1], tmp1$group))[2, ] / sum(tmp1[, 1])
    pct <- rbind.data.frame(pct, tmp1)
  }
  rownames(pct) <- lesion.order
  
  # Right-side stacked percentage bar plot
  right_anno2 <- anno_barplot(as.matrix(pct),
                              which = "row",
                              border = FALSE,
                              gp = gpar(fill = clust.col,
                                        border = NA,
                                        lty = "blank"), 
                              bar_width = 0.6,
                              width = unit(1.8, "cm"),
                              height = unit(1, "cm"))
  # Draw the heatmap in the same way
  op2 <- oncoPrint(onco.input2[lesion.order, rownames(my_ann)], 
                   alter_fun = alter_fun2, 
                   col = col2, 
                   bottom_annotation = NULL, 
                   top_annotation = my_annotation, # Top annotation: subtype
                   column_order = rownames(my_ann),
                   right_annotation = rowAnnotation(PCT = right_anno2),
                   row_order = lesion.order, 
                   show_pct = T,
                   column_title = "", 
                   show_heatmap_legend = T, 
                   column_split = my_ann$group,
                   column_title_gp = gpar(fontsize = 8),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   alter_fun_is_vectorized = FALSE)
  op2
  pdf("mutational arm level in TCGA.pdf", width = 8, height = 6)
  draw(op2) 
  invisible(dev.off())
  ### Start chi-square test heatmap
  tmp_aow <- tmp %>% filter(group != "High")
  tmp_High <- tmp %>% filter(group == "High")
  ## Try
  matrix(c(sum(tmp_aow[,1] == "1"), sum(tmp_aow[,1] == "0"), sum(tmp_High[,1] == "1"), sum(tmp_High[,1] == "0")), nrow = 2)
  chisq.test(matrix(c(sum(tmp_aow[,1] == "1"), sum(tmp_aow[,1] == "0"), sum(tmp_High[,1] == "1"), sum(tmp_High[,1] == "0")), nrow = 2)
  )
  ## Batch calculate CNV frequency and rate
  chiq_arm <- data.frame(arm = colnames(tmp)[1:(ncol(tmp)-1)])
  for (i in 1:(ncol(tmp)-1)) {
    summary = chisq.test(matrix(c(sum(tmp_aow[,chiq_arm$arm[i]] == "1"), sum(tmp_aow[,chiq_arm$arm[i]] == "0"), sum(tmp_High[,chiq_arm$arm[i]] == "1"), sum(tmp_High[,chiq_arm$arm[i]] == "0")), nrow = 2)
    )
    p = summary$p.value
    High = sum(tmp_High[,chiq_arm$arm[i]] == "1")
    Low = sum(tmp_aow[,chiq_arm$arm[i]] == "1")
    chiq_arm$p_value[i] = p
    chiq_arm$`High group`[i] = High
    chiq_arm$`Low group`[i] = Low
    chiq_arm$High[i] = round(High / length(rownames(tmp_High)), 2)
    chiq_arm$Low[i] = round(Low / length(rownames(tmp_aow)), 2)
  }
  ## Start plotting
  B <- chiq_arm[, c("arm", "High", "Low")] %>% pivot_longer(cols = !arm, 
                                                            names_to = "Cohort",
                                                            values_to = "frequency")
  cnvarm <- character()
  for (i in 1:length(lesion.order)) { 
    g = lesion.order[i]
    if (chiq_arm$p_value[chiq_arm == g] < 0.001) {
      cnvarm[i] = paste0("*** ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] < 0.01) {
      cnvarm[i] = paste0("** ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] < 0.05) {
      cnvarm[i] = paste0("* ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] > 0.05) {
      cnvarm[i] = g
    }
  }
  for (i in 1:length(B$arm)) {
    g = B$arm[i]
    if (chiq_arm$p_value[chiq_arm == g] < 0.001) {
      B$arm[i] = paste0("*** ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] < 0.01) {
      B$arm[i] = paste0("** ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] < 0.05) {
      B$arm[i] = paste0("* ", g)
    } else if (chiq_arm$p_value[chiq_arm == g] > 0.05) {
      B$arm[i] = g
    }
  }
  B$arm <- factor(B$arm, levels = rev(cnvarm))
  
  B$Cohort <- factor(B$Cohort, levels = c("Low", "High"))
  p1 <- ggplot(B, aes(x = Cohort, y = arm, fill = `frequency`)) +
    
    geom_raster() +
    geom_tile(lwd = 1,
              linetype = 1, colour = "black") +
    theme_minimal() +
    scale_fill_gradient(limits = c(0, 1), low = 'white', high = 'red') +
    geom_text(aes(label = sprintf("%.3f", `frequency`)), col = 'black', cex = 4) +
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title = element_blank()) + 
    scale_y_discrete(position = "left")
  
  p1
  
  group_heatmap <- ggplot(data = B) + geom_tile(aes(x = Cohort, y = 1, fill = Cohort), lwd = 1,
                                                linetype = 1, colour = "black") +
    scale_fill_manual(values = c(Low = '#87CEEB', High = '#EE2C2C')) +
    
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), axis.line = element_blank(), 
          axis.ticks = element_blank(), axis.text = element_blank(), 
          axis.title = element_blank(),
          legend.position = "right") 
  group_heatmap
  p_arm <- p1 %>% insert_top(group_heatmap, height = .05)
  print(p_arm)
  export::graph2pdf(file = "heatmap_arm.pdf", width = 3.7, height = 5)
  dev.off()
  
  ### Batch plot arm
  gene_matrix <- test_cna[lesion.order, rownames(my_ann)] %>% rbind(target, .) %>% t() %>% as.data.frame()  
  gene_matrix = gather(gene_matrix, "cnvarm", "Group", 2:ncol(gene_matrix))
  gene_matrix$Group = factor(gene_matrix$Group, c("Not CNV", "CNV"))
  gene_matrix$cnvarm = factor(gene_matrix$cnvarm, `lesion.order`) ## Fix gene order
  gene_matrix$gene <- as.numeric(gene_matrix$gene)
  colnames(gene_matrix)[1] = paste0(genedds, " Expression")
  p <- ggboxplot(gene_matrix, x = "Group", y = paste0(genedds, " Expression"), color = "Group", 
                 palette = "jco", # Use jco palette, can also change to "lancet"
                 add = "jitter", facet.by = "cnvarm", short.panel.labs = T,
                 size = 0.3, width = 0.5) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, label.y = 3) # R defaults to method = "wilcox.test" # Add p-value, run this line to get the figure below
  
  print(p)
  export::graph2pdf(file = paste0(genedds, paste0(paste0(" TCGA-", ca), "batch cnv arms.pdf")), width = 8, height = 8)
  dev.off()
  
  # rm(list = setdiff(ls(), c("try", "cas")))
}
for (ca in cas) {try(ca)
  
}
rm(list = ls())
