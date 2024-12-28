library(Seurat)
library(ggplot2)

#---load object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
obj$sample <- factor(obj$sample, levels = rev(c("wt-12w", "dare-12w")))

#---print violin plot---
a <- VlnPlot(obj, features = "Pdgfra", assay = "RNA", pt.size = 0, group.by = "annotation",
             split.by = "sample") 
a + scale_fill_manual(
  values = c(
    "wt-12w" = "grey70",
    "dare-12w" = "red2"
  ),
  labels = c(
    "wt-12w" = expression(italic('Tnf'^'+/+')),
    "dare-12w" = expression(italic('Tnf'^'DARE'))
  )
) + theme(
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 45), 
  legend.text.align = 0,
  text = element_text(size = 22)
)
ggsave(filename = "SFig12e.pdf", width = 25, height = 10, dpi = 600, units = "cm")
