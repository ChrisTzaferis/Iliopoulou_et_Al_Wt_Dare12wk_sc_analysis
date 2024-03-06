library(Seurat)

#Load objects
dareLym <- readRDS("../Step5_DoubletFinder/dare12_Lymphoid_PCs_23_0.05_DoubletFinder.RDS")
dareMye <- readRDS("../Step5_DoubletFinder/dare12_Myeloid_PCs_21_0.5_DoubletFinder.RDS")
dareStr <- readRDS("../Step5_DoubletFinder/dare12_Stroma_PCs_22_0.05_DoubletFinder.RDS")

wtLym <- readRDS("../Step5_DoubletFinder/wt12_Lymphoid_PCs_15_0.05_DoubletFinder.RDS")
wtMye <- readRDS("../Step5_DoubletFinder/wt12_Myeloid_PCs_14_0.5_DoubletFinder.RDS")
wtStr <- readRDS("../Step5_DoubletFinder/wt12_Stroma_PCs_20_0.05_DoubletFinder.RDS")

integrateSeuratObjects <- function(obj1, obj2, nPCs, resolution, exportName)
{
  #change Idents and keep only unique cells
  Idents(obj1) <- obj1$DoubletFinderClass
  Idents(obj2) <- obj2$DoubletFinderClass
  obj1 <- subset(obj1, idents = "FALSE")
  obj2 <- subset(obj2, idents = "FALSE")
  
  #pre-processing 
  obj.list <- list(obj1, obj2)
  # normalize and identify variable features for each dataset independently
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = obj.list)
  
  #Perform integration
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
  
  # this command creates an 'integrated' data assay
  obj.combined <- IntegrateData(anchorset = obj.anchors)
  
  #Perform integrated analysis
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(obj.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  obj.combined <- ScaleData(obj.combined, verbose = T)
  obj.combined <- RunPCA(obj.combined, npcs = 30, verbose = T)
  ElbowPlot(obj.combined, ndims = 30)
  obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:nPCs)
  obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:nPCs)
  obj.combined <- FindClusters(obj.combined, resolution = 0.5)
  
  #save RDS
  saveRDS(object = obj.combined, file = exportName)
}

integrateSeuratObjects(obj1 = wtStr, obj2 = dareStr, nPCs = 23, resolution = 0.5, exportName = "WT12wk_DARE12Wk_Integration_Stroma_VST_without_DoubletsFromDoubletFinder.RDS")
integrateSeuratObjects(obj1 = wtLym, obj2 = dareLym, nPCs = 15, resolution = 0.5, exportName = "WT12wk_DARE12Wk_Integration_Lymphoid_VST_without_DoubletsFromDoubletFinder.RDS")
integrateSeuratObjects(obj1 = wtMye, obj2 = dareMye, nPCs = 20, resolution = 0.5, exportName = "WT12wk_DARE12Wk_Integration_Myeloid_VST_without_DoubletsFromDoubletFinder.RDS")
