library(Seurat)

wt12 <- readRDS("../Step3_SeuratWorkflowPerSample/WT12_weeks_before_integration_nF_800.RDS")
dare12 <- readRDS("../Step3_SeuratWorkflowPerSample/DARE12_weeks_before_integration_nF_800.RDS")

analyzeAndSaveObjects <- function(seu_obj, supercluster, nPCs, resolution, sampleName)
{
  obj_filtered <- subset(seu_obj, idents = supercluster)
  
  #normalization
  obj_filtered <- NormalizeData(obj_filtered)
  
  #HVGs
  obj_filtered <- FindVariableFeatures(obj_filtered, selection.method = "mvp")
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(obj_filtered), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(obj_filtered)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  #scaling
  #all.genes <- rownames(obj_filtered)
  obj_filtered <- ScaleData(obj_filtered, vars.to.regress = c("nCount_RNA", "percent.mt"))
  
  #PCA
  obj_filtered <- RunPCA(obj_filtered, features = VariableFeatures(object = obj_filtered))
  ElbowPlot(obj_filtered, ndims = 50)
  
  #Clustering
  obj_filtered <- FindNeighbors(obj_filtered, dims = 1:nPCs)
  obj_filtered <- FindClusters(obj_filtered, resolution = resolution)
  
  #UMAP
  obj_filtered <- RunUMAP(obj_filtered, dims = 1:nPCs)
  obj_filtered <- RunTSNE(obj_filtered, dims = 1:nPCs)
  DimPlot(obj_filtered, reduction = "umap")
  
  #export RDS
  saveRDS(object = obj_filtered, file = paste0(sampleName, "_", supercluster, "_PCs_", nPCs, "_", resolution, ".RDS")) 
}

analyzeAndSaveObjects(seu_obj = wt12, supercluster = "Myeloid", nPCs = 14, resolution = 0.5, sampleName = "wt12")
analyzeAndSaveObjects(seu_obj = dare12, supercluster = "Myeloid", nPCs = 21, resolution = 0.5, sampleName = "dare12")

analyzeAndSaveObjects(seu_obj = wt12, supercluster = "Lymphoid", nPCs = 15, resolution = 0.5, sampleName = "wt12")
analyzeAndSaveObjects(seu_obj = dare12, supercluster = "Lymphoid", nPCs = 23, resolution = 0.5, sampleName = "dare12")

analyzeAndSaveObjects(seu_obj = wt12, supercluster = "Stroma", nPCs = 20, resolution = 0.5, sampleName = "wt12")
analyzeAndSaveObjects(seu_obj = dare12, supercluster = "Stroma", nPCs = 22, resolution = 0.5, sampleName = "dare12")
