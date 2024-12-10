library(Seurat)
library(ggplot2)

#---load object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

#---print violin plot---
a <- VlnPlot(obj, features = "Cd14", assay = "RNA", pt.size = 0.1, idents = c("Resident_Mfs", "Activated_Mfs"), 
        split.by = "sample") 
a + scale_fill_manual(
          values = c(
            "wt-12w" = "grey70",
            "dare-12w" = "red"
          ),
          labels = c(
            "wt-12w" = expression(italic('Tnf'^'+/+')),
            "dare-12w" = expression(italic('Tnf'^'DARE'))
          )
        ) + theme(
                  axis.title.x = element_blank(), 
                  axis.text.x = element_text(angle = 0, hjust = 0.5), 
                  legend.text.align = 0,
                  text = element_text(size = 22)
                  )
ggsave(filename = "Fig2e.pdf", width = 16, height = 10, dpi = 600, units = "cm")
