library(Seurat)
library(Scillus)
library(ggplot2)

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
  scale_x_discrete(labels=c("wt-12w" = expression(italic("Tnf"^"+/+")), "dare-12w" = expression(italic("Tnf"^"DARE")) )) +
  theme_classic() + 
  coord_flip() +
  theme(legend.position = "bottom", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 18)
        )
ggsave("SFig3b.pdf", device = "pdf", width = 24, height = 7, units = "cm", dpi = 600)
