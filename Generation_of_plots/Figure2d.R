library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

#---load myeloid object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                 "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

#---Dotplot---
genes <- c("Cxcl2", "Il23a", "Ccl3", "Ccl4", "Il1b", "Il1a", 
           "Il1rn", "Tnf", "Cxcl1", "Ccl5", "Il6", "Osm", "Il12a", "Lif")
DotPlot(obj, features = genes, dot.scale = 13, group.by = "annotation", assay = "RNA", scale = T) +
  theme_classic() +
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  theme(text = element_text(size = 25), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="grey")))
ggsave(filename = "Fig2d.pdf", device = "pdf", width = 26, height = 11, units = "cm", dpi = 600)
