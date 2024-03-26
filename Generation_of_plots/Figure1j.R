library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(corrplot)
library(viridis)
library(pheatmap)

int_obj <- readRDS("../Step8_Preparation_of_human_dataset/lymphoid_integration.RDS")
#--keep dare and inflamed samples
Idents(int_obj) <- int_obj$sample
int_obj <- subset(int_obj, idents = c("dare-12w", "Infl"))

#--correlation input
DefaultAssay(object = int_obj) <- "integrated"

int_obj$cluster_sample <- paste0(int_obj$sample, "_", int_obj$annotation)
Idents(int_obj) <- int_obj@meta.data$cluster_sample

avgExp <- AverageExpression(object = int_obj, assays = "integrated", slot = "data")
avgExp_Array <- as.data.frame(avgExp$integrated)
avgExp_Array$Gene <- rownames(avgExp_Array)

gene_col <- grep('Gene', x = colnames(avgExp_Array))
avgDataMat <- as.matrix(avgExp_Array[, -gene_col] )
rownames(avgDataMat) <- avgExp_Array$Gene

#--calculate correlation
res <- cor(avgDataMat, method = "spearman")
res <- round(res, 2)

#--row and numbers in the final matrix
all_clusters <- length(colnames(res)) 
dare_clusters <- length(grep(pattern = "^dare", x = colnames(res)))
infl_clusters <- all_clusters - dare_clusters

#--matrix initiation
sp_corr_mat <- matrix(0, nrow = infl_clusters, ncol = dare_clusters)
colnames(sp_corr_mat) <- colnames(res)[grep(pattern = "^dare", x = colnames(res))]
rownames(sp_corr_mat) <- colnames(res)[grep(pattern = "^Infl", x = colnames(res))]

#--fill matrix
for(i in 1:nrow(sp_corr_mat))
  for(j in 1:ncol(sp_corr_mat))
  {
    sp_corr_mat[i, j] <- res[which(rownames(res) == rownames(sp_corr_mat)[i]), which(colnames(res) == colnames(sp_corr_mat)[j])]
  }

#--print heatmap
pdf("Fig1j.pdf", width = 10, height = 10)
pheatmap(
  sp_corr_mat[, ],
  color = colorRampPalette(c("blue", "white", "red"))(82),
  breaks = seq(0.52, 0.92, by=0.005),
  cluster_cols = T,
  cluster_rows = T, display_numbers = F,
  fontsize = 10
)
dev.off()
