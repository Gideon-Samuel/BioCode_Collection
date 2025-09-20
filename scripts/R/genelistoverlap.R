# R Script: Venn Diagram for Gene Lists
# - Visualize overlap between two gene sets

#Installing and loading packages
if(!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)

#Example gene lists
# In practice, these could come from DE analysis results of two conditions
genes_A <- c("TP53","MYC","BRCA1","EGFR","PTEN","CDK2")
genes_B <- c("BRCA1","EGFR","PTEN","MTOR","CDK4","KRAS")

#Finding overlaps
common_genes <- intersect(genes_A, genes_B)
cat("Common genes:", common_genes, "\n")

#Venn Diagram
venn.plot <- draw.pairwise.venn(area1 = length(genes_A),
                                area2 = length(genes_B),
                                cross.area = length(common_genes),
                                category = c("List A", "List B"),
                                fill = c("lightblue", "pink"),
                                alpha = 0.5,
                                cat.pos = c(-20, 20),
                                cat.dist = 0.05)

#Save Venn as PNG

png("gene_overlap.png", width=600, height=600)
draw.pairwise.venn(area1 = length(genes_A),
                   area2 = length(genes_B),
                   cross.area = length(common_genes),
                   category = c("List A", "List B"),
                   fill = c("lightblue", "pink"),
                   alpha = 0.5,
                   cat.pos = c(-20, 20),
                   cat.dist = 0.05)
dev.off()
