# R Script: Airway RNA-seq Dataset Exploration
# - Human airway smooth muscle cells ± dexamethasone
# - Visualize expression of a selected gene

# Installing required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("airway", quietly = TRUE)) BiocManager::install("airway")

# Loading libraries
library(airway)
library(ggplot2)
library(dplyr)

# Loading the dataset
data("airway")

# Extract the expression matrix and metadata
air_data <- as.data.frame(assay(airway))
meta <- as.data.frame(colData(airway))       

# Checking sample names
colnames(air_data)
meta$dex

gene_id <- "ENSG00000103196"

# Building the expression dataframe for plotting
gene_expr <- data.frame(
  Sample = colnames(airway_data),
  Condition = meta$dex,
  Expression = as.numeric(airway_data[gene_id, ])
)

# Plotting
ggplot(gene_expr, aes(x = Condition, y = Expression, color = Condition)) +
  geom_jitter(width = 0.1, size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  ylab(paste0(gene_id, " Expression")) +
  ggtitle(paste("Expression of", gene_id, "across samples"))
