# R Script: Explore GEO Dataset GSE47472 (Microarray)
# - Human samples (details in GEO)
# - QC: boxplot, density, PCA for sample clustering

# Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")

# Loading libraries
library(GEOquery)
library(limma)
library(affy)
library(ggplot2)

#Downloading GEO dataset
gse <- getGEO("GSE47472", GSEMatrix = TRUE) 
expr <- exprs(gse[[1]])
meta <- pData(gse[[1]])

#Boxplot for sample distributions
boxplot(expr,
        main = "Boxplot of Microarray Expression (GSE47472)",
        col = "lightblue", las = 2,
        outline = FALSE)

#Density plot for expression distribution per sample)
matplot(density(expr[,1])$x,
        sapply(1:ncol(expr), function(i) density(expr[,i])$y),
        type = "l", lty = 1, col = rainbow(ncol(expr)),
        xlab = "Expression", ylab = "Density",
        main = "Density Plot of Expression (QC)")
legend("topright", legend = colnames(expr), col = rainbow(ncol(expr)), lty = 1, cex = 0.5)

#PCA plot for sample clustering
pca <- prcomp(t(expr), scale. = TRUE)   # PCA on samples
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     Sample = rownames(pca$x),
                     Group = meta$characteristics_ch1.1) # may differ by dataset

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  ggtitle("PCA of Microarray Samples (GSE47472)")
