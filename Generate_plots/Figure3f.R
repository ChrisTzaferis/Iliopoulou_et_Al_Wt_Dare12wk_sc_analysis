source("utilityFunctions/correlationAnalysisSeuratObject.R")

obj <- readRDS("../Step8_Preparation_of_human_dataset/stroma_integration.RDS") #seurat object

correlationHeatmap(obj, usedAssay = "integrated", usedCondition = "sample", usedSamples = c("Infl", "dare-12w"), 
                   usedClusters = "annotation", usedSlot = "data", corMethod = "spearman", usedGenes = VariableFeatures(obj), 
                   objFileName = "Fig3f", clustCols = F, clustRows = F, 
                   row_order = c(2, 5, 8, 11, 14, 9, 16, 12, 3, 7, 10, 15, 6, 4, 13, 1), 
                   col_order = c(9, 14, 5, 8, 7, 11, 2, 12, 3, 6, 1, 10, 4, 13)
                   )
