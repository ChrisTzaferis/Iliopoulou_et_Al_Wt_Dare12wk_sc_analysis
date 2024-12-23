library(Scillus)
library(ggplot2)
library(Seurat)

#load object and order clusters, samples
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
obj$seurat_clusters <- obj$annotation
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

cols <- c("skyblue1","slateblue4", "salmon1", "yellow4", "thistle", "yellow3", 
                       "red2","yellowgreen", "mediumorchid2",
                       "wheat3", "royalblue", "orange","grey70","tan3")

#Barplot
plot_stat(obj, plot_type = "prop_fill", group_by = "sample", pal_setup = cols) + 
  theme_classic() + 
  scale_x_discrete(labels=c("wt-12w" = expression(italic("Tnf"^"+/+")), "dare-12w" = expression(italic("Tnf"^"DARE")) )) +
  theme(legend.position = "right", 
        axis.title.x = element_blank(),
        text = element_text(size = 18)
  ) + 
  labs(y = "Percent of cells")
ggsave("SFig5c.pdf", device = "pdf", width = 15, height = 25, units = "cm", dpi = 600)
