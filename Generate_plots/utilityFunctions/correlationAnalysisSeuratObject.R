library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(corrplot)
library(viridis)
library(pheatmap)

#Parameters
#obj <- readRDS("...RDS") #a seurat object
#usedAssay <- "integrated" #RNA or integrated
#usedCondition <- "sample" #column containing conditions or NULL for correlating just clusters
#usedSamples <- c("dare-12w", "inflamed") #pair of samples if usedCondition is not NULL
#usedClusters <- "annotation" #column containing active clusters in metadata
#usedSlot <- "data" #counts, data, scale.data
#corMethod <- "pearson" #"pearson", "kendall", "spearman"
#usedGenes <- NULL #keep a subset of genes, NULL by default -> all genes are used
#objFileName <- "Fig1j" #text in the filename of the exported pdf
#clustRows #clustering of Rows (TRUE/FALSE)
#clustCols #clustering of Cols (TRUE/FALSE)
#row_order #custom ordering of rows, (Print heatmap without clustering the first time and then modify order)
#col_order #custom ordering of columns, (Print heatmap without clustering the first time and then modify order)

correlationHeatmap <- function(obj, usedAssay="RNA", usedCondition=NULL, usedSamples=c("", ""), usedClusters="seurat_sclusters", usedSlot="data", 
                               corMethod = "spearman", usedGenes=NULL, objFileName="", clustRows=T, clustCols=T, 
                               row_order=NULL, col_order=NULL)
{
  #--subset selected samples if needed
  if(!is.null(usedCondition) )
    if(usedCondition %in% colnames(obj@meta.data))
    {
      Idents(obj) <- usedCondition
      obj <- subset(obj, idents = usedSamples)
    }
  
  #--correlation input
  DefaultAssay(obj) <- usedAssay
  
  #
  if(!is.null(usedCondition) )
  {
    if(usedCondition %in% colnames(obj@meta.data))
    {
      obj$cluster_sample <- paste0(obj@meta.data[[usedClusters]], "_", obj@meta.data[[usedCondition]])
      Idents(obj) <- obj@meta.data$cluster_sample
    }  
  }else
  {
    Idents(obj) <- usedClusters
  }
  
  avgExp <- AverageExpression(object = obj, assays = usedAssay, slot = usedSlot)
  avgExp_Array <- as.data.frame(avgExp[[usedAssay]])
  avgExp_Array$Gene <- rownames(avgExp_Array)
  
  #keep only genes of interest or all genes
  if(!is.null(usedGenes))
  {
    avgExp_Array <- avgExp_Array[which(avgExp_Array$Gene %in% usedGenes), ]
  }
  
  #convert to matrix
  gene_col <- grep('Gene', x = colnames(avgExp_Array))
  avgDataMat <- as.matrix(avgExp_Array[, -gene_col] )
  rownames(avgDataMat) <- avgExp_Array$Gene
  
  #--calculate correlation
  res <- cor(avgDataMat, method = corMethod)
  res <- round(res, 2)
  
  #--row and numbers in the final matrix
  all_clusters <- length(colnames(res)) 
  
  if(!is.null(usedCondition))
  {
    if(usedCondition %in% colnames(obj@meta.data))
    {
      cond1_clusters <- length(grep(pattern = paste0("_", usedSamples[1], "$"), x = colnames(res)))
      cond2_clusters <- all_clusters - cond1_clusters
      
      #--matrix initiation
      sp_corr_mat <- matrix(0, nrow = cond1_clusters, ncol = cond2_clusters)
      rownames(sp_corr_mat) <- colnames(res)[grep(pattern = paste0("_", usedSamples[1], "$"), x = colnames(res))]
      colnames(sp_corr_mat) <- colnames(res)[grep(pattern = paste0("_", usedSamples[2], "$"), x = colnames(res))]
      
      #--fill matrix
      for(i in 1:nrow(sp_corr_mat))
        for(j in 1:ncol(sp_corr_mat))
        {
          sp_corr_mat[i, j] <- res[which(rownames(res) == rownames(sp_corr_mat)[i]), which(colnames(res) == colnames(sp_corr_mat)[j])]
        }
    }
  }else
  {
    sp_corr_mat <- res
  }
  
  #printing order in rows and cols of heatmap
  col_ordering <- 1:ncol(sp_corr_mat)
  if(clustCols == F & !is.null(col_order))
  {
    col_ordering <- col_order
  }
  
  row_ordering <- 1:nrow(sp_corr_mat)
  if(clustRows == F & !is.null(row_order))
  {
    row_ordering <- row_order
  }
  
  #--heatmap colors
  # Choose color palette
  my_breaks <- seq(min(sp_corr_mat), max(sp_corr_mat), by=0.005)
  palette_length = length(my_breaks)+1
  my_color = colorRampPalette(c("blue", "white","red"))(palette_length)
  
  #--print heatmap
  filename <- paste0(objFileName, ".pdf")
  pdf(filename, width = 10, height = 10)
  pheatmap(
    sp_corr_mat[row_ordering, col_ordering],
    color = my_color,
    breaks = my_breaks,
    cluster_cols = clustCols,
    cluster_rows = clustRows, 
    display_numbers = F,
    fontsize = 10
  )
  dev.off()
}
