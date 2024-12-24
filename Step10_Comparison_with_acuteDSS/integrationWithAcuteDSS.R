library(Seurat)
library(scCustomize)
library(dplyr)
library(dittoSeq)
library(tidyverse)

#acute DSS obj
obj_dss <- readRDS("DSSacute_HC_integrated_annotated.RDS")
obj_dss$sample <- obj_dss$orig.ident
obj_dss$annotation <- Idents(obj_dss)

#dare obj
obj_dare <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")

DimPlot(obj_dss)
DimPlot(obj_dare)

#First keep only stroma cells for both datasets
obj_dare <- subset(obj_dare, idents = c("Pdgfra_low", "Telocytes", "Trophocytes", "FRC1", "FRC2"))
DimPlot(obj_dare, label = T, split.by = "sample")

obj_dss <- subset(obj_dss, idents = c("S1", "S2", "S3", "Frc"))
DimPlot(obj_dss, label = T, split.by = "sample")

DimPlot(obj_dare)
DimPlot(obj_dss)

#view metadata
tail(obj_dare@meta.data)
tail(obj_dss@meta.data)

#keep the same metadata columns
obj_dare@meta.data <- obj_dare@meta.data[, c(1, 2, 3, 5, 21)]
obj_dss@meta.data <- obj_dss@meta.data[, c(1, 2, 3, 6, 7)]


#Integration analysis
DefaultAssay(obj_dss) <- "RNA"
obj_dss[['integrated']] <- NULL
DefaultAssay(obj_dare) <- "RNA"
obj_dare[['integrated']] <- NULL

dare_list <- SplitObject(obj_dare, split.by = "sample")
dss_list <- SplitObject(obj_dss, split.by = "sample")

obj.list <- list(dare_list[[1]], dss_list[[1]], dare_list[[2]], dss_list[[2]])
# identify variable features for each dataset independently
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures=2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.list)

#Perform integration
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# this command creates an 'integrated' data assay
obj.combined <- IntegrateData(anchorset = obj.anchors)

#Perform integrated analysis on the correct assay
DefaultAssay(obj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.combined <- ScaleData(obj.combined, verbose = T)
obj.combined <- RunPCA(obj.combined, npcs = 30, verbose = T)
ElbowPlot(obj.combined, ndims = 30)
nPCs <- 20
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:nPCs)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:nPCs)
obj.combined <- FindClusters(obj.combined, resolution = 0.5)

DimPlot_scCustom(obj.combined, split.by = "sample", group.by = "annotation", label = T, pt.size = 0.8, label.size = 6,
                 colors_use = scCustomize_Palette(num_groups = 9, ggplot_default_colors = F)[-2])

#Rename Frcs and save object
obj.combined$annotation <- gsub(x = obj.combined$annotation, pattern = "Frc", replacement = "S4")
obj.combined$sample <- factor(obj.combined$sample, levels = c("wt-12w", "HC", "dare-12w", "DSS_acute"))

saveRDS(obj.combined, "wt12_HC_dare12_DSSacute_Stroma_integration.RDS")