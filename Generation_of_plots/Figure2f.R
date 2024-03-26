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
cols <- c("cornflowerblue", "aquamarine4", "darkgoldenrod", "red3")
Idents(obj) <- obj$sample
wt <- subset(obj, idents = "wt-12w")
dare <- subset(obj, idents = "dare-12w")

#Wt umap
p1 <- dittoDimPlot(wt, var = "annotation", color.panel = cols, size = 1.4, show.axes.numbers = F, show.grid.lines = F, 
             add.trajectory.lineages = list(
               c("Monocytes","Intermediate_Mfs", "Resident_Mfs"),
               c("Monocytes","Intermediate_Mfs","Activated_Mfs")
             ), trajectory.cluster.meta = "annotation"
             ) + ggtitle("") + theme_void() + theme(legend.position = "none") + xlim(-1, 10) + ylim(-7, 4)

p2 <- dittoDimPlot(dare, var = "annotation", color.panel = cols, size = 1.4, show.axes.numbers = F, show.grid.lines = F, 
             add.trajectory.lineages = list(
               c("Monocytes","Intermediate_Mfs", "Resident_Mfs"),
               c("Monocytes","Intermediate_Mfs","Activated_Mfs")
             ), trajectory.cluster.meta = "annotation"
) + ggtitle("") + theme_void() + theme(legend.position = "none") + xlim(-1, 10) + ylim(-7, 4)

cowplot::plot_grid(p1, p2,
                   labels = c('Wt','DARE'),
                   ncol = 1, align = "v", label_size = 24)

ggsave("Fig2f.pdf", device = "pdf", width = 10, height = 25, units = "cm", dpi = 600)

#Barplot
obj$seurat_clusters <- obj$annotation
obj$seurat_clusters <- factor(obj$seurat_clusters, levels = c("Activated_Mfs", "Intermediate_Mfs", "Monocytes", "Resident_Mfs"))
plot_stat(obj, plot_type = "prop_fill", group_by = "sample", pal_setup = c("red3", "aquamarine4", "cornflowerblue", "darkgoldenrod")) + 
  theme_classic() + 
  theme(legend.position = "right", 
        text = element_text(size = 16),
        axis.title.x = element_blank()) + 
  ylab("Percent of cells") +
  ggtitle("")
ggsave("Fig2f_barplot.pdf", device = "pdf", width = 15, height = 25, units = "cm", dpi = 600)
