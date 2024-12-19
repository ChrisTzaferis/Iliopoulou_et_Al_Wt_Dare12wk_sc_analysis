library(Scillus)
library(ggplot2)
library(Seurat)

#load object and order clusters, samples
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")
final_order<-c("B_cells","Cycling_B_cells","Plasma_cells", "Naïve_CD4_T_cells", "Naïve_CD8_T_cells","Cd8_T_cells", "Memory_T_cells",
               "Th17_cells","Tregs", "Tcr_gammadelta", "Cycling_T_cells","NKT_cells","ILC2s","ILC3s")
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$seurat_clusters <- obj$annotation
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

cols <- c("#CEB2DF", "#E4A149", "#7B78D9", "#B940E1", "#D9A99C", "#D5E2D6", 
          "#70B3D6", "#D76DC6", "#DC637B", "#70E16E", "#7BE9CE", "#D0E24B",
          "#D2DE95", "#7C8E66")

#Barplot
plot_stat(obj, plot_type = "prop_fill", group_by = "sample", pal_setup = cols) + 
  theme_classic() + 
  scale_x_discrete(labels=c("wt-12w" = expression(italic("Tnf"^"+/+")), "dare-12w" = expression(italic("Tnf"^"DARE")) )) +
  theme(legend.position = "right", 
        axis.title.x = element_blank(),
        text = element_text(size = 18)
        ) + 
  labs(y = "Percent of cells")
ggsave("SFig2b.pdf", device = "pdf", width = 15, height = 25, units = "cm", dpi = 600)
