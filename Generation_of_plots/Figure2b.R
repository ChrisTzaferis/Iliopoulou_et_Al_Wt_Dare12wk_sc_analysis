library(Seurat)
library(scCustomize)
library(ggplot2)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
DefaultAssay(obj) <- "RNA"

#Feature plots
FeaturePlot_scCustom(obj, features = c("Adgre1", "Cx3cr1", "Zbtb46", "Flt3", "S100a9", "S100a8"), num_columns = 2,
                     na_cutoff = NA, colors_use = c("grey70", "red"), pt.size = 1) & 
  theme(legend.position = "right", line = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave("Fig2b.pdf", device = "pdf", width = 17, height = 22, units = "cm", dpi = 600)
