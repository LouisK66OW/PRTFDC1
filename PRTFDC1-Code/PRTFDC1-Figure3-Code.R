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

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################### KM Analysis

load(file = "KM.Rdata")
#surv<-TCGA_PFI
#surv<-OS_3218
#surv<-OS_10783
# Determine the optimal cutoff value for gene expression
cutoff <- surv_cutpoint(surv, # Dataset
                        time="time", # Time variable in the dataset
                        event="status", # Outcome variable name in the dataset
                        variables="ENSG00000099256",
                        minprop = 0.3 # Column name to be calculated
);summary(cutoff) # Output results

# Visualize the cutoff value of a gene's expression level
plot(cutoff, 
     "ENSG00000099256", 
     palette = "lancet") # Use Lancet color palette
groups <- surv_categorize(cutoff)
str(groups)
head(groups)
fit <- survfit(Surv(time, status) ~ ENSG00000099256, data=groups) # Survival analysis function

# Plot survival curve
ggsurvplot(fit,
           data = surv,
           conf.int = F,
           pval.method=TRUE,  # Show p-value calculation method
           risk.table = TRUE,
           surv.median.line = "hv",
           xlab="Time(months)",  # Set x-axis name
           ylab="Overall Survival",  # Set y-axis name
           break.time.by = 50,  # Set x-axis interval
           title="", # Main title
           legend.labs=c("High-risk","Low-risk"),  # Add legend labels
           legend.title="PRTFDC1", # Add legend title
           palette = "hue",
           font.legend = c(20, "plain"),      # Legend font settings
           font.x = c(20, "plain", "black"),  # x-axis label font settings
           font.y = c(20, "plain", "black"),  # y-axis label font settings
           axes.offset = TRUE, # Whether axes origin intersect
           legend = c(0.9, 0.9), # Change legend position
           pval = T, pval.size = 8, pval.coord = c(100, .25), # p-value font size and position
           pval.method.coord = c(100, .35),
           font.xtickslab = c(20, "plain", "black"),  # x-axis tick label font settings
           font.ytickslab = c(20, "plain", "black"),  # y-axis tick label font settings
           risk.table.title = "Number at risk", # risk.table title
           risk.table.fontsize = 10, # risk.table font size
           risk.table.col = "strata", # risk.table text color, consistent with the curve
           risk.table.y.text = TRUE, # Display text on y-axis of risk.table, FALSE for color blocks
           tables.y.text = FALSE,
           risk.table.height = 0.25, # Proportion of risk.table in the image
           surv.plot.height = 0.75, # Proportion of survival plot in the image
           font.main = c(22, "bold")
) # Main title font settings

rm(list = ls())

load(file = "COX.Rdata")

TGCT_m<-cox_TCGA_PFI
#TGCT_m<-cox_OS_3218
#TGCT_m<-cox_OS_10783

TGCT_m <- cox_TCGA_PFI
## Batch univariate COX analysis
## Adjust data to avoid errors (convert gene_symbol)
data = TGCT_m
genes <- colnames(data)[-c(1:2)]
library(survival)
res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv = as.formula(paste('Surv(time, status) ~', genes[i]))
  y = coxph(surv, data = data)
  x = summary(y)
  p.value = signif(x$wald["pvalue"], digits=3)
  HR = signif(x$coef[2], digits=3); # exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 3)
  HR.confint.upper = signif(x$conf.int[,"upper .95"], 3)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
  res[i,5] = HR.confint.lower
  res[i,6] = HR.confint.upper
}

names(res) <- c("Variable","HR","95% CI","p.value","HR.confint.lower","HR.confint.upper")

res0.05 <- res[res$p.value < 0.05,]
res0.2 <- res[res$p.value < 0.2,]
res0.5 <- res[res$p.value < 0.5,]

library(survival)
### Multivariate COX regression analysis
multicox <- coxph(Surv(time = time, event = status) ~ ., data = TGCT_m) 
multisum <- summary(multicox)
# Extract results of multivariate COX regression analysis for all genes into multiresult object
gene <- colnames(TGCT_m)[3:ncol(TGCT_m)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(gene=gene,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)

library(survival)
library(survminer)
ggforest(multicox,
         data = TGCT_m,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 1.5,
         refLabel = "reference",
         noDigits = 3
)

