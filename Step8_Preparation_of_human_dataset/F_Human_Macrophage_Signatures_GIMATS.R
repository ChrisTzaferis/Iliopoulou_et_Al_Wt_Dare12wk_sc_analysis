library(scCustomize)
library(Seurat)
library(dplyr)
library(stringr)
library(nichenetr)

#Files were downloaded from GEO [Series GSE134809] and only ileum samples were utilized for further analysis 
all_dirs <- list.dirs("GIMATS_files")
all_dirs <- all_dirs[-c(1, 22)]
all_dirs

#Create a merged Seurat object
my_list <- list()

for(i in 1:length(all_dirs))
{
  current_sample <- all_dirs[i]
  proj_name_list <- str_split(current_sample, "/", n = 2, simplify = T)
  proj_name <- proj_name_list[2]
  
  temp <- Read10X(current_sample)
  temp_obj <- CreateSeuratObject(counts = temp, project = proj_name, min.features = 200)
  temp_obj$sample <- proj_name
  my_list <- append(my_list, temp_obj)
}

merged <- merge(my_list[[1]], c(my_list[[2]], my_list[[3]], my_list[[4]], my_list[[5]], my_list[[6]], my_list[[7]], my_list[[8]], my_list[[9]], my_list[[10]], 
                                my_list[[11]], my_list[[12]], my_list[[13]], my_list[[14]], my_list[[15]], my_list[[16]], my_list[[17]], my_list[[18]], 
                                my_list[[19]], my_list[[20]]), merge.data = F)
merged
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged@meta.data$merged_annot <- "MergedIleum"
saveRDS(merged, "mergedIleum_unfiltered.RDS")

#-------------------------------------------------------------------------------
#Filtering
epithelial_list <- c("PLA2G2A", "CLCA1", "REG4", "S100A14", "ITLN1", "ELF3", "PIGR", "EPCAM", "REG1B", "REG1A", "REG3A", "FABP1", "RBP2", 
                     "SST", "FABP2", "SPINK1", "FABP6", "AGR2", "AGR3", "CLDN3", "CLDN4", "DEFA6", "DEFA5", "SPINK4", "ALDOB", "LCN2", "MUC2", 
                     "KRT8", "KRT18", "TSPAN8", "OLFM4", "GPX2", "IFI27", "PHGR1", "MT1G", "CLDN7", "KRT19", "FXYD3", "LGALS4", "FCGBP", 
                     "TFF3", "TFF1")
merged[["percent.epi"]] <- PercentageFeatureSet(merged, features = epithelial_list[which(epithelial_list %in% rownames(merged))])

#merged #108734 cells
obj <- subset(merged, subset = percent.epi <= 1 & percent.mt <= 25 & nCount_RNA >= 800) #72872 cells (using filtering described in "Methods" section of the publication)
obj
saveRDS(obj, "GIMATS_ileum_filtered.RDS")

#-------------------------------------------------------------------------------
#load and store metadata info 
obj <- readRDS("GIMATS_ileum_filtered.RDS")
metadata <- read.delim("metadataGIMATS.txt")
metadata$sample <- paste0(metadata$tissue, "_", metadata$Sample_ID)
metadata$sample <- gsub(x = metadata$sample, pattern = "ILEUM", replacement = "Ileal")
obj@meta.data$Cell_id <- rownames(obj@meta.data) 
obj@meta.data <- left_join(obj@meta.data, metadata)
rownames(obj@meta.data) <- obj$Cell_id

#Remove patient 6
Idents(obj) <- obj$Patient.ID
obj <- subset(obj, idents = "pat. 6", invert=T)
table(obj$Patient.ID)

#Default analysis
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "merged_annot")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "merged_annot")
plot1 + plot2

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.epi"), ncol = 4, group.by = "merged_annot", pt.size = 0)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

obj <- ScaleData(obj)

obj <- RunPCA(obj, features = VariableFeatures(object = obj))

ElbowPlot(obj, ndims = 50)

obj <- FindNeighbors(obj, dims = 1:14)
obj <- FindClusters(obj, resolution = 0.5)

obj <- RunUMAP(obj, dims = 1:14)

DimPlot(obj, reduction = "umap", label = T, label.size = 5, split.by = "status")#, group.by = "orig.ident")

FeaturePlot_scCustom(obj, features = c("ITGAX", "ITGAM", "PTPRC"))
FeaturePlot_scCustom(obj, features = c("ITGAX"), colors_use = viridis_light_high, pt.size = 0.9)
FeaturePlot_scCustom(obj, features = c("HLA-DRB1"), colors_use = viridis_light_high, pt.size = 0.9)
FeaturePlot_scCustom(obj, features = c("LYZ"), colors_use = viridis_light_high, pt.size = 0.9)
saveRDS(obj, "analyzedGIMATS_allCells.RDS")
#-------------------------------------------------------------------------------
# Extract myeloid cells
obj <- readRDS("analyzedGIMATS_allCells.RDS")
DimPlot_scCustom(obj, group.by = "RNA_snn_res.0.5", label = T)
obj_sub <- subset(obj, idents = c("6", "14", "17"))
DimPlot(obj_sub)

#Re-analysis
obj <- obj_sub
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

obj <- ScaleData(obj)

obj <- RunPCA(obj, features = VariableFeatures(object = obj))

ElbowPlot(obj, ndims = 50)

obj <- FindNeighbors(obj, dims = 1:25)
obj <- FindClusters(obj, resolution = c(0.1, 0.2))

obj <- RunUMAP(obj, dims = 1:25)

DimPlot_scCustom(obj, reduction = "umap", label = T, label.size = 5, split.by = "status", group.by = "RNA_snn_res.0.2")
saveRDS(obj, "myeloid_subclusters.RDS")

################################################################################
#annotate clusters based on marker genes
# 0 Resident_Macs
# 1 DC2
# 2 Inflammatory_Macs
# 3 Activated_DC
# 4 MoDC
# 5 DC1
# 6 Resident_Macs
# 7 pDC
# 8 DC2

Idents(obj) <- obj$RNA_snn_res.0.2
obj <- Rename_Clusters(seurat_object = obj, new_idents = c("Resident_Macs", "DC2", "Inflammatory_Macs", "Activated_DC", 
                                                           "MoDC", "DC1", "Resident_Macs", "pDC", "DC2"))
obj$annotation <- Idents(obj)
DimPlot(obj, label = T, group.by = c("annotation", "RNA_snn_res.0.2"), 
        cols = scCustomize_Palette(num_groups = 9, ggplot_default_colors = F))
saveRDS(obj, "myeloid_subclusters_annotated.RDS")