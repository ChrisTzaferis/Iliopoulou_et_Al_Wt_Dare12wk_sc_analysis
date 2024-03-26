library(Seurat)
library(scCustomize)
library(ggplot2)

#load objects
wt <- readRDS("../Step3_SeuratWorkflowPerSample/WT12_weeks_before_integration_nF_800.RDS")
dare <- readRDS("../Step3_SeuratWorkflowPerSample/DARE12_weeks_before_integration_nF_800.RDS")

#UMAP plots
DimPlot(wt, cols = c("purple3", "coral", "grey70"), label = T, pt.size = 1, label.size = 10) + 
  theme_classic() + 
  theme(legend.position = "none", line = element_blank(), text = element_blank(), title = element_blank())
ggsave("Fig1b_umap_wt.pdf", device = "pdf", width = 20, height = 20, units = "cm", dpi = 600)


DimPlot(dare, cols = c("purple3", "grey70", "coral"), label = T, pt.size = 1, label.size = 10) + 
  theme_classic() + 
  theme(legend.position = "none", line = element_blank(), text = element_blank(), title = element_blank())
ggsave("Fig1b_umap_dare.pdf", device = "pdf", width = 20, height = 20, units = "cm", dpi = 600)

#Feature plots
FeaturePlot_scCustom(wt, features = c("Itgax", "Itgam", "Ptprc"), num_columns = 1,
                     na_cutoff = NA, colors_use = c("grey70", "red"), pt.size = 1) & 
  theme(legend.position = "none", line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), 
        panel.border = element_rect(color = "grey20"))
ggsave("Fig1b_FPs_normal.pdf", device = "pdf", width = 12, height = 22, units = "cm", dpi = 600)


FeaturePlot_scCustom(dare, features = c("Itgax", "Itgam", "Ptprc"), num_columns = 1,
                     na_cutoff = NA, colors_use = c("grey70", "red"), pt.size = 1) & 
  theme(legend.position = "none", line = element_blank(), axis.text = element_blank(), axis.title = element_blank(), 
        panel.border = element_rect(color = "grey20"))
ggsave("Fig1b_FPs_dare.pdf", device = "pdf", width = 12, height = 22, units = "cm", dpi = 600)