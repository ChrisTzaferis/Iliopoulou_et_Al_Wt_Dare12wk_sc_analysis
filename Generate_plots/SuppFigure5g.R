library(ggplot2)
library(pheatmap)
library(Seurat)
library(dplyr)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")

#keep clusters of interest
obj <- subset(obj, idents = c("Pdgfra_low", "Trophocytes", "Telocytes"))

#change clustering
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)
Idents(obj) <- obj$cluster_sample

#selected genes
gene_list <- data.frame(Gene = c("Col6a1", "Col6a2", "Col6a3", "Fn1", "Col1a2", 
               "Col1a1", "Col3a1", "Col2a1", "Col5a1", "Col15a1"))

#calculate average expression
avg_table <- AverageExpression(obj, assays = "RNA")
avg_table <- avg_table$RNA
avg_df <- as.data.frame(avg_table)
avg_df$Gene <- rownames(avg_df)
avg_htbale <- left_join(gene_list, avg_df)
avg_htbale <- na.omit(avg_htbale)
avg_hmat <- as.matrix(avg_htbale[, -1])
rownames(avg_hmat) <- avg_htbale$Gene

#pheatmap
mat_scaled <- t(scale(t(avg_hmat+1)))

pdf("SFig5g.pdf", width = 20, height = 25)
pheatmap(mat_scaled,
         show_rownames = T,
         cluster_cols = T, 
         cellwidth = 120, 
         cellheight = 120,
         fontsize = 34, 
         legend_breaks = c(-1, 0, 1),
         color = colorRampPalette(c("blue", "white", "red"))(50)
         )
dev.off()
