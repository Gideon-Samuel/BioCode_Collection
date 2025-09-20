# R Script: Summary Statistics and Visualization for GSE70970 RNA-seq - Human dataset
# - Compute per-gene statistics and visualize expression
# Installing required packages (Bioconductor/CRAN)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

# Loading libraries
library(GEOquery)
library(ggplot2)
library(dplyr)
library(tidyr)

#Loading packages
library(GEOquery)
library(ggplot2)
library(dplyr)

#Downloading GEO dataset
gse <- getGEO("GSE70970", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
meta <- pData(gse[[1]])

#Computing summary statistics per gene
gene_stats <- data.frame(
  Gene = rownames(expr),
  Mean = rowMeans(expr),
  SD = apply(expr, 1, sd),
  Median = apply(expr, 1, median),
  IQR = apply(expr, 1, IQR)
)

# Viewing first few genes
head(gene_stats)

#Histogram of mean expression
ggplot(gene_stats, aes(x = Mean)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
  xlab("Mean Expression") + ylab("Number of Genes") +
  ggtitle("Distribution of Mean Gene Expression")

# Boxplot of expression for first 20 genes
expr_subset <- expr[1:20, ]
expr_long <- as.data.frame(t(expr_subset))
expr_long$Sample <- rownames(expr_long)
expr_long <- tidyr::pivot_longer(expr_long, cols = -Sample, names_to = "Gene", values_to = "Expression")

ggplot(expr_long, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Expression of First 20 Genes Across Samples")

