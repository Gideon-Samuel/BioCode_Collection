# R Script: Differential Expression Analysis (GSE19804)
# - Lung cancer samples: Stage 1A vs Stage 3A
# - Limma for DE, volcano plot for visualization


#Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

#Loading libraries
library(GEOquery)
library(limma)
library(ggplot2)

#Downloading GEO dataset
gse <- getGEO("GSE19804", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
meta <- pData(gse[[1]])
head(meta)

#Defining group: stage 1A vs 3A
meta$Group <- ifelse(meta$`stage:ch1` == "1A", "Stage1A",
                     ifelse(meta$`stage:ch1` == "3A", "Stage3A", NA))

# Removing samples not in 1A or 3A
keep <- !is.na(meta$Group)
expr <- expr[, keep]
meta <- meta[keep, ]
group <- factor(meta$Group)
table(group)

#Designing matrix for DE

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#Fit linear model and computing DE

fit <- lmFit(expr, design)
contrast <- makeContrasts(Stage3A - Stage1A, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
de_results <- topTable(fit2, number = nrow(expr), adjust.method = "BH")

#Volcano plotting

de_results$Significant <- "No"
de_results$Significant[de_results$adj.P.Val < 0.05 & abs(de_results$logFC) > 1] <- "Yes"

ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot: Stage3A vs Stage1A (GSE19804 Lung Cancer)") +
  theme_minimal()
