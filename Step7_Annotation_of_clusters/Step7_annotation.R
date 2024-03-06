library(Seurat)
library(scCustomize)
library(ggplot2)

#-----------------------Annotation of STROMA cells------------------------------
stroma <- readRDS("../Step6_Integration/WT12wk_DARE12Wk_Integration_Stroma_VST_without_DoubletsFromDoubletFinder.RDS")
DefaultAssay(stroma) <- "integrated"

#Exclusion of clusters 10, 14, 19
stroma <- subset(x = stroma, idents = c(10, 14, 19), invert = T)
stroma <- FindNeighbors(stroma, reduction = "pca", dims = 1:23)
stroma <- FindClusters(stroma, resolution = 0.4)

#Annotation of clusters after manual inspection of marker genes from the literature
new.cluster.ids <- c("Pdgfra_low","Pdgfra_low", "LECS", "Trophocytes", "Telocytes", "SMCs","FRC2", "FRC1", "Capillary_ECs", "Trophocytes", "Pericytes", 
                     "FDCs_MRCs","Cajal", "Vein_ECs", "Mesothelial", "Glial")
names(new.cluster.ids) <- levels(stroma)
stroma <- RenameIdents(object = stroma, new.cluster.ids)
stroma$annotation <- Idents(stroma)
saveRDS(stroma, "WT12wk_DARE12Wk_Stroma_annotated.RDS")
#-------------------------------------------------------------------------------


#-----------------------Annotation of Lymphoid cells----------------------------
lymphoid <- readRDS("../Step6_Integration/WT12wk_DARE12Wk_Integration_Lymphoid_VST_without_DoubletsFromDoubletFinder.RDS")
DefaultAssay(lymphoid) <- "integrated"
lymphoid <- FindClusters(lymphoid, resolution = 0.7)
DimPlot_scCustom(lymphoid)

#Exclusion of cluster 16
lymphoid <- subset(lymphoid, idents = 16, invert = T)
DimPlot_scCustom(lymphoid)
lymphoid
lymphoid <- FindNeighbors(lymphoid, reduction = "pca", dims = 1:23)
lymphoid <- FindClusters(lymphoid, resolution = 0.7)

#Annotation of clusters after manual inspection of marker genes from the literature
new.cluster.ids <- c("B_cells","Plasma_cells","Plasma_cells","Plasma_cells", "Memory_T_cells", "ILC2s", "Naïve_CD4_T_cells", 
                     "Cd8_T_cells", "Th17_cells", "ILC3s", "Tregs", "Cycling_B_cells", "Tcr_gammadelta", "B_cells","Naïve_CD8_T_cells",
                     "Cycling_T_cells", "Plasma_cells","NKT_cells","Memory_T_cells", "B_cells")
names(new.cluster.ids) <- levels(lymphoid)
lymphoid <- RenameIdents(lymphoid, new.cluster.ids)
lymphoid$annotation <- Idents(lymphoid)
saveRDS(lymphoid, "WT12wk_DARE12Wk_Lymphoid_annotated.RDS")
#-------------------------------------------------------------------------------


#-----------------------Annotation of Myeloid cells-----------------------------
myeloid <- readRDS("../Step6_Integration/WT12wk_DARE12Wk_Integration_Myeloid_VST_without_DoubletsFromDoubletFinder.RDS")
DefaultAssay(myeloid) <- "integrated"
myeloid <- FindClusters(myeloid, resolution = 0.8)
DimPlot_scCustom(myeloid)

#Exclusion of clusters 11, 12, 13, 14, 16
myeloid <- subset(myeloid, idents = c("11", "12", "13", "14", "16"), invert = T)
myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:23)
myeloid <- FindClusters(myeloid, resolution = 0.5)

#Annotation of clusters after manual inspection of marker genes from the literature
new.cluster.ids <- c("Intermediate_Mfs","Granulocytes", "Granulocytes", "Activated_Mfs", "Resident_Mfs",
                     "Cdc2_a","Monocytes", "Cdc1", "Mixed_Mf_popul", "Cdc2_b", "pDCs", "Cdc2_b")
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)

#exclusion of a mixed MF population
myeloid <- subset(myeloid, idents = c("Mixed_Mf_popul"), invert = T)
DimPlot_scCustom(myeloid)
myeloid$annotation <- Idents(myeloid)
saveRDS(myeloid, "WT12wk_DARE12Wk_Myeloid_annotated.RDS")
#-------------------------------------------------------------------------------