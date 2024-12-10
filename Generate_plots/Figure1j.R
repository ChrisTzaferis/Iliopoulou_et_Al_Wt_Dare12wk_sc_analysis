source("./utilityFunctions/correlationAnalysisSeuratObject.R")

obj <- readRDS("../Step8_Preparation_of_human_dataset/lymphoid_integration.RDS") #a seurat object
correlationHeatmap(obj, usedAssay = "integrated", usedCondition = "sample", usedSamples = c("Infl", "dare-12w"), 
                   usedClusters = "annotation", usedSlot = "data", corMethod = "spearman", usedGenes = VariableFeatures(obj), 
                   objFileName = "Fig1j", clustRows = T, clustCols = T)
