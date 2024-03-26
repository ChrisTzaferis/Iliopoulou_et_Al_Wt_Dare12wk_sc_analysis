library(Seurat)
library(Scillus)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Activated_Mfs", "Cdc1", "Cdc2_a", "Cdc2_b", "Granulocytes", "Intermediate_Mfs", "Monocytes", "pDCs", "Resident_Mfs")
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$seurat_clusters <- obj$annotation
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

cols <- c("red3","darkgreen",
          "hotpink3", "lightpink3", "chocolate2", 
          "aquamarine4", "cornflowerblue","ivory4", "darkgoldenrod")

#Barplot
plot_stat(obj, plot_type = "prop_fill", group_by = "sample", pal_setup = cols) + 
  theme_classic() + 
  coord_flip() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
 ggtitle("Percent of cells")
ggsave("Fig2d_barplot.pdf", device = "pdf", width = 20, height = 5, units = "cm", dpi = 600)
