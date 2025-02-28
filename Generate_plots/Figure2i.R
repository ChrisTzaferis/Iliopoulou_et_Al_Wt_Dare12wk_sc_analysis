library(Seurat)
library(dittoSeq)
library(Scillus)

#Load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                 "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

#keep macrophages 
obj <- subset(obj, idents = c("Monocytes", "Intermediate_Mfs", "Activated_Mfs", "Resident_Mfs"))
Idents(obj) <- obj$sample
wt <- subset(obj, idents = "wt-12w")
dare <- subset(obj, idents = "dare-12w")
cols <- c("#D4D8DA", "#B2DFEE", "#DDCC77", "#882255")

#Wt umap
p1 <- dittoDimPlot(wt, var = "annotation", color.panel = cols, size = 1.4, show.axes.numbers = F, show.grid.lines = F, 
             add.trajectory.lineages = list(
               c("Monocytes","Intermediate_Mfs", "Resident_Mfs"),
               c("Monocytes","Intermediate_Mfs","Activated_Mfs")
             ), trajectory.cluster.meta = "annotation"
             ) + ggtitle(expression(italic("Tnf"^"+/+"))) + 
                 theme_void() + 
                 theme(legend.position = "none", text = element_text(size = 22)) + 
                 xlim(-1, 10) + ylim(-7, 4)

p2 <- dittoDimPlot(dare, var = "annotation", color.panel = cols, size = 1.4, show.axes.numbers = F, show.grid.lines = F, 
             add.trajectory.lineages = list(
               c("Monocytes","Intermediate_Mfs", "Resident_Mfs"),
               c("Monocytes","Intermediate_Mfs","Activated_Mfs")
             ), trajectory.cluster.meta = "annotation"
             ) + ggtitle(expression(italic("Tnf"^"DARE"))) + 
                 theme_void() + 
                 theme(legend.position = "none", text = element_text(size = 22)) + 
                 xlim(-1, 10) + ylim(-7, 4)

cowplot::plot_grid(p1, p2,
                   ncol = 1, align = "v", label_size = 24)

ggsave("Fig2i_umap.pdf", device = "pdf", width = 10, height = 25, units = "cm", dpi = 600)

#Barplot
obj$seurat_clusters <- obj$annotation
obj$seurat_clusters <- factor(obj$seurat_clusters, levels = c("Activated_Mfs", "Resident_Mfs", "Intermediate_Mfs", "Monocytes"))
plot_stat(obj, plot_type = "prop_fill", group_by = "sample", pal_setup = rev(cols)) + 
  scale_x_discrete(labels=c("wt-12w" = expression(italic("Tnf"^"+/+")), "dare-12w" = expression(italic("Tnf"^"DARE")) )) +
  theme_classic() + 
  theme(legend.position = "right", 
        text = element_text(size = 18),
        axis.title.x = element_blank()) + 
  ylab("Percent of cells") +
  ggtitle("")
ggsave("Fig2i_barplot.pdf", device = "pdf", width = 15, height = 25, units = "cm", dpi = 600)
