library(Seurat)
library(scCustomize)
library(dplyr)

#GSE114374 download data from GEO
#Load DSS, HC matrix
mouse_hc_mat <- read.delim("GSE114374_Mouse_HC_expression_matrix.txt")
mouse_dss_mat <- read.delim("GSE114374_Mouse_DSS_expression_matrix.txt")

#Create Seurat objects
obj1 <- CreateSeuratObject(counts = mouse_hc_mat, project = "HC")
#max(obj1@assays$RNA@data)

obj2 <- CreateSeuratObject(counts = mouse_dss_mat, project = "DSS_acute")
#max(obj2@assays$RNA@data)

# Integration analysis
obj1
obj2

obj.list <- list(obj1, obj2)
# identify variable features for each dataset independently (counts were already normalized)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.list)

#Perform integration
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# this command creates an 'integrated' data assay
obj.combined <- IntegrateData(anchorset = obj.anchors)

#Perform integrated analysis in the correct assay
DefaultAssay(obj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.combined <- ScaleData(obj.combined, verbose = T)
obj.combined <- RunPCA(obj.combined, npcs = 30, verbose = T)
ElbowPlot(obj.combined, ndims = 30)
nPCs <- 20
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:nPCs)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:nPCs)
obj.combined <- FindClusters(obj.combined, resolution = 0.3)

DimPlot_scCustom(obj.combined, split.by = "orig.ident", group.by = "integrated_snn_res.0.3", label = T, pt.size = 0.8, label.size = 6,
                 colors_use = scCustomize_Palette(num_groups = 14, ggplot_default_colors = F)[-2])

#Check gene expression patterns
DefaultAssay(obj.combined) <- "RNA"
Idents(obj.combined) <- obj.combined$integrated_snn_res.0.3
FeaturePlot(obj.combined, features = c("Ccl19", "Il11", "Il33"), ncol = 3, order = T, label = T)

FeaturePlot(obj.combined, features = c("Myh11", "Ano1", "Pecam1", "Rgs5", "Lyve1", "Cdh19"), 
            ncol = 3, order = T, label = T)

FeaturePlot(obj.combined, features = c("Ptprc", "Epcam"), pt.size = 0.8,
            ncol = 2, order = T, label = F)

#S1
FeaturePlot(obj.combined, features = c("Pdgfra", "Negr1", "Vcan", "Edil3"), 
            ncol = 3, order = T, label = T)

#S2
FeaturePlot(obj.combined, features = c("Bmp7", "Bmp5", "Sox6", "Foxl1"), 
            ncol = 3, order = T, label = T)

#S3
FeaturePlot(obj.combined, features = c("Col14a1", "Gpc3", "Slit2", "Cd81"), 
            ncol = 3, order = T, label = T)

#S4
FeaturePlot(obj.combined, features = c("Ccl19", "Il11", "Il33", "Ccl19"), 
            ncol = 3, order = T, label = T)

#Rename clusters
Idents(obj.combined)
obj.combined <- Rename_Clusters(seurat_object = obj.combined, new_idents = c("S3", "S1", "S1", "S2", "S3", 
                                                                             "Frc", "Smcs", "Lecs", "EndCs", "pericytes", 
                                                                             "Frc", "Immune_cells", "Glia"))
DimPlot_scCustom(obj.combined, split.by = "orig.ident")

saveRDS(obj.combined, "DSSacute_HC_integrated_annotated.RDS")
