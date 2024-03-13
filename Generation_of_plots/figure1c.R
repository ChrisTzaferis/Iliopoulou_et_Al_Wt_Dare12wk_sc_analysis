library(Seurat)
library(scCustomize)
library(ggplot2)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")
final_order<-c("B_cells","Cycling_B_cells","Plasma_cells", "Naïve_CD4_T_cells", "Naïve_CD8_T_cells","Cd8_T_cells", "Memory_T_cells",
               "Th17_cells","Tregs", "Tcr_gammadelta", "Cycling_T_cells","NKT_cells","ILC2s","ILC3s")
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

#UMAP plot
DimPlot_scCustom(obj, group.by = "annotation", colors_use = c("#CEB2DF", "#E4A149", "#7B78D9", "#B940E1", "#D9A99C", "#D5E2D6", 
                                                              "#70B3D6", "#D76DC6", "#DC637B", "#70E16E", "#7BE9CE", "#D0E24B",
                                                              "#D2DE95", "#7C8E66"
), 
                 split.by = "sample", 
                 figure_plot = T, 
                 split_seurat = T,
                 label = T, 
                 pt.size = 0.8, 
                 repel = T
)
ggsave("Fig1c.pdf", device = "pdf", width = 30, height = 20, units = "cm", dpi = 600)

