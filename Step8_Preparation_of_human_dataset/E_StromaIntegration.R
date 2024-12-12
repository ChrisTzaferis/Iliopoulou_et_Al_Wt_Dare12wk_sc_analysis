library(Seurat)

#load object from immunity paper
obj1 <- readRDS("TI_Stroma_Human_with_mouse_genes.RDS")
head(obj1@meta.data)
obj1$Cell_id <- obj1$NAME
obj1$sample <- obj1$Type
obj1$annotation <- obj1$Celltype
head(obj1@meta.data)
Idents(obj1) <- obj1$sample

#load DARE+WT object
obj2 <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
head(obj2@meta.data)
obj2$sample <- factor(obj2$sample, levels = c("wt-12w", "dare-12w"))
Idents(obj2) <- obj2$sample
DefaultAssay(obj2) <- "RNA"

##########################--------------------##################################
#########################    Integration     ###################################

obj.list <- list(subset(obj1, idents = "Infl"),
                 subset(obj1, idents = "NonI"),
                 subset(obj1, idents = "Heal"),
                 subset(obj2, idents = "wt-12w"),
                 subset(obj2, idents = "dare-12w")
)

# normalize and identify variable features for each dataset independently
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.list)

obj.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)

# this command creates an 'integrated' data assay
obj.combined <- IntegrateData(anchorset = obj.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(obj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.combined <- ScaleData(obj.combined, verbose = T)
obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = T)
ElbowPlot(obj.combined, ndims = 50)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:22)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:22)
obj.combined <- FindClusters(obj.combined, resolution = 1.2)

# Visualization
DimPlot(obj.combined, reduction = "umap", group.by = "annotation", raster = F)
saveRDS(obj.combined, "stroma_integration.RDS")