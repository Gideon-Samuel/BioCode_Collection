# R Script: Explore GEO Dataset GSE18388 (mouse thymus, spaceflight study)
# Expression matrix and sample metadata extraction


# Installing GEOquery
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")

#Loading the Library
library(GEOquery)

# Downloading the GEO dataset
gse <- getGEO("GSE18388", GSEMatrix = TRUE)

# Extracting expression data
expr <- exprs(gse[[1]])

# Extracting phenotype data
meta <- pData(gse[[1]])

# Viewing at dataset
dim(expr)          
head(expr[, 1:5])
head(meta)         
View(meta)