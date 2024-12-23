library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

#Load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
final_order<-c("Pdgfra_low","LECS","Trophocytes", "Telocytes", "SMCs","FRC2", "FRC1",
               "Capillary_ECs","Pericytes", "FDCs_MRCs", "Cajal","Vein_ECs","Mesothelial","Glial")
obj$annotation <- factor(obj$annotation, levels = final_order)
DefaultAssay(obj) <- "RNA"

#---Markers---
all_markers <- FindAllMarkers(obj, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T, base = exp(1))

#---Dotplot---
cluster_markers <-  all_markers %>% arrange(desc(avg_logFC)) %>%
  group_by(cluster) %>% 
  slice(1:10)

DotPlot(obj, features = unique(cluster_markers$gene), dot.scale = 13, group.by = "annotation", assay = "RNA", scale = T) +
  theme_classic() +
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  theme(text = element_text(size = 25), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )+
  coord_flip()+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="grey")))
ggsave(filename = "SFig5b.pdf", device = "pdf", width = 30, height = 115, units = "cm", dpi = 600)
