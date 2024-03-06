library(Seurat)
library(patchwork)
library(dplyr)

#object1 initialization
obj1_counts <- Read10X("../GEO_files/PROCESSED/WT12/")
obj1 <- CreateSeuratObject(counts = obj1_counts, project = "WT12wk", min.cells = 0, min.features = 0)
obj1[["percent.mt"]] <- PercentageFeatureSet(obj1, pattern = "^mt-")
obj1$sample <- "wt-12w"
obj1$genotype <- "wt"

#doublet removal
#load scrublet results
obj1_doublets <- read.csv("Step1_DoubletRemovalWithScrublet/WT12_scrublet_results.csv")
colnames(obj1_doublets)[1] <- "Cell_id"
obj1_doublets <- obj1_doublets[, c("Cell_id", "predicted_doublets")]

metadata_temp1 <- as.data.frame(obj1@meta.data)
metadata_temp1$Cell_id <- rownames(metadata_temp1)
metadata_temp1 <- left_join(metadata_temp1, obj1_doublets)
rownames(metadata_temp1) <- metadata_temp1$Cell_id

obj1@meta.data <- metadata_temp1
Idents(obj1) <- obj1$predicted_doublets

#View QC plots before doublet removal
VlnPlot(obj1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0)

#Remove doublets
obj1_no_doublets <- subset(obj1, idents = "False")

obj1 #43119 cells
obj1_no_doublets #35887 cells

#View QC plots after doublet removal
VlnPlot(obj1_no_doublets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0)
saveRDS(obj1_no_doublets, "wt12_no_doubletsDefaultThreshold_not_filtered.RDS")

#object2 initialization
obj2_counts <- Read10X("../GEO_files/PROCESSED/DARE12/")
obj2 <- CreateSeuratObject(counts = obj2_counts, project = "DARE12wk", min.cells = 0, min.features = 0)
obj2[["percent.mt"]] <- PercentageFeatureSet(obj2, pattern = "^mt-")
obj2$sample <- "dare-12w"
obj2$genotype <- "dare"

#doublet removal
#load scrublet results
obj2_doublets <- read.csv("Step1_DoubletRemovalWithScrublet/DARE12_scrublet_results.csv")
colnames(obj2_doublets)[1] <- "Cell_id"
obj2_doublets <- obj2_doublets[, c("Cell_id", "predicted_doublets")]

metadata_temp2 <- as.data.frame(obj2@meta.data)
metadata_temp2$Cell_id <- rownames(metadata_temp2)
metadata_temp2 <- left_join(metadata_temp2, obj2_doublets)
rownames(metadata_temp2) <- metadata_temp2$Cell_id

obj2@meta.data <- metadata_temp2
Idents(obj2) <- obj2$predicted_doublets

#View QC plots before doublet removal
VlnPlot(obj2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0)

#Remove doublets
obj2_no_doublets <- subset(obj2, idents = "False")

obj2 #25771 cells
obj2_no_doublets #22708 cells

#View QC plots after doublet removal
VlnPlot(obj2_no_doublets, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident", pt.size = 0)
saveRDS(obj2_no_doublets, "dare12_no_doubletsDefaultThreshold_not_filtered.RDS")
