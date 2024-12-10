library(Seurat)
library(scCustomize)
library(ggplot2)
library(patchwork)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")
final_order<-c("B_cells","Cycling_B_cells","Plasma_cells", "Naïve_CD4_T_cells", "Naïve_CD8_T_cells","Cd8_T_cells", "Memory_T_cells",
               "Th17_cells","Tregs", "Tcr_gammadelta", "Cycling_T_cells","NKT_cells","ILC2s","ILC3s")
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

obj_list <- SplitObject(obj, split.by = "sample")

#UMAP plot for WT
g1 <- DimPlot_scCustom(obj_list[[1]], group.by = "annotation", colors_use = c("#CEB2DF", "#E4A149", "#7B78D9", "#B940E1", "#D9A99C", "#D5E2D6", 
                                                              "#70B3D6", "#D76DC6", "#DC637B", "#70E16E", "#7BE9CE", "#D0E24B",
                                                              "#D2DE95", "#7C8E66"), 
                 figure_plot = F, 
                 split_seurat = F,
                 label = T, 
                 pt.size = 0.8, 
                 repel = T
) + theme(legend.position = "none", 
          line = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25)) + 
  labs(title = expression(italic('Tnf'^'+/+')))


#UMAP plot for DARE
g2 <- DimPlot_scCustom(obj_list[[2]], group.by = "annotation", colors_use = c("#CEB2DF", "#E4A149", "#7B78D9", "#B940E1", "#D9A99C", "#D5E2D6", 
                                                                              "#70B3D6", "#D76DC6", "#DC637B", "#70E16E", "#7BE9CE", "#D0E24B",
                                                                              "#D2DE95", "#7C8E66"), 
                       figure_plot = F, 
                       split_seurat = F,
                       label = T, 
                       pt.size = 0.8, 
                       repel = T
) + theme(legend.position = "right", 
          line = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25)) + 
  labs(title = expression(italic('Tnf'^'DARE')))

#Combination of plots with patchwork
g1 + g2 + plot_annotation(title = 'Myeloid', theme = theme(plot.title = element_text(size = 20, hjust = 0.4)))
ggsave("Fig1c.pdf", device = "pdf", width = 30, height = 20, units = "cm", dpi = 600)
