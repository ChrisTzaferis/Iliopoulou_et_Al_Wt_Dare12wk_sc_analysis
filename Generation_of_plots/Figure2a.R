library(Seurat)
library(scCustomize)
library(ggplot2)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                 "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

#UMAP plot
DimPlot_scCustom(obj, group.by = "annotation", colors_use = c("chocolate2","cornflowerblue",
                                                              "aquamarine4", "darkgoldenrod", "red3", 
                                                              "darkgreen", "hotpink3","lightpink3", "ivory4"
), 
split.by = "sample", 
figure_plot = T, 
split_seurat = T,
label = T, 
pt.size = 0.8, 
repel = T
)
ggsave("Fig2a.pdf", device = "pdf", width = 30, height = 20, units = "cm", dpi = 600)
