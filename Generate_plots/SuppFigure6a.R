library(Seurat)
library(ggplot2)
library(scCustomize)
library(tidyverse)

#load integrated object
obj.combined <- readRDS("../Step10_Comparison_with_acuteDSS/wt12_HC_dare12_DSSacute_Stroma_integration.RDS")

#umap plots
obj.combined$annotation <- factor(obj.combined$annotation)
g2 <- DimPlot_scCustom(obj.combined, 
                 pt.size = 0.8, 
                 reduction = "umap", 
                 group.by = "annotation", 
                 split.by = "sample",  
                 split_seurat = F,
                 label = F,
                 colors_use = c("yellowgreen", "mediumorchid2", "tan3", "chocolate4", "darkgoldenrod3", 
                                "grey40", "forestgreen", "orange", "grey70")
)

obj.combined@meta.data <- obj.combined@meta.data |> mutate(dataset = ifelse(sample %in% c("HC", "DSS_acute"), "DSS_dataset", "DARE_dataset") )
g1 <- DimPlot_scCustom(obj.combined, 
                 pt.size = 0.8, 
                 reduction = "umap", 
                 group.by = "dataset", 
                 split_seurat = T, 
                 shuffle = T,
                 label = F,
                 colors_use = c("cyan4", "burlywood")
)

#Modify plots before save
integ_plot <- g1 + 
  theme(legend.position = "bottom") +
  labs(title = "") +
  scale_color_manual(
    values = c(
      "DSS_dataset" = "burlywood",
      "DARE_dataset" = "cyan4"
    ),
    labels = c(
      "DSS_dataset" = expression(bold("DSS_dataset")), 
      "DARE_dataset" = expression(bolditalic("Tnf"^{DARE})~bold("dataset"))
    )) 

wt_plot <- g2[[1]] + ggtitle(label = "Healthy ileum")
hc_plot <- g2[[2]] + ggtitle(label = "Healthy colon")
dare_plot <- g2[[3]] + ggtitle(label = expression(bolditalic("Tnf"^{DARE})~bold("ileum")) )
dss_plot <- g2[[4]] + ggtitle(label = "DSS colon")
integ_plot + ((wt_plot + hc_plot) / (dare_plot + dss_plot))
ggsave(filename = "SFig6a.pdf", width = 45, height = 20, device = "pdf", dpi = 600, units = "cm")
