library(Seurat)
library(ggplot2)
library(scCustomize)

#load object
obj <- readRDS("../Step2_InitialSeuratObjects/wt12_no_doubletsDefaultThreshold_not_filtered.RDS")
obj

#QC filtering
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
obj_filtered <- subset(obj, subset = nFeature_RNA > 800 & nFeature_RNA < 7500 & percent.mt < 5)
VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
obj_filtered

#wt: 
# before filtering: 35428 cells
# nFeature_RNA > 800 & nFeature_RNA < 7500 & percent.mt < 5: 24249 cells

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
obj_filtered <- FindNeighbors(obj_filtered, dims = 1:20)
obj_filtered <- FindClusters(obj_filtered, resolution = 0.01)

#UMAP
obj_filtered <- RunUMAP(obj_filtered, dims = 1:20)
DimPlot(obj_filtered, reduction = "umap", label = T)

#annotation of WT
new.cluster.ids <- c("Lymphoid", "Stroma", "Lymphoid", "Myeloid", "Stroma", "Stroma")
names(new.cluster.ids) <- levels(obj_filtered)
obj_filtered <- RenameIdents(obj_filtered, new.cluster.ids)

#save
saveRDS(obj_filtered, "WT12_weeks_before_integration_nF_800.RDS") 

#used in Fig1b
FeaturePlot(obj_filtered, features = "Itgax", order = T, label = T, cols = c("grey70", "red"))
FeaturePlot(obj_filtered, features = "Itgam", order = T, label = T, cols = c("grey70", "red"))
FeaturePlot(obj_filtered, features = "Ptprc", order = T, label = T, cols = c("grey70", "red"))
FeaturePlot(obj_filtered, features = "Col1a1", order = T, label = T, cols = c("grey70", "red"))

DimPlot_scCustom(obj_filtered, colors_use = c("purple", "coral", "grey"), pt.size = 0.8, label.size = 7) + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank()
        )
