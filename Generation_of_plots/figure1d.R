library(Seurat)
library(ggplot2)
library(scDataviz)

#load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")

#Contour plot
Idents(obj) <- obj$sample
DefaultAssay(obj) <- "RNA"
sce1 <- as.SingleCellExperiment(subset(obj, idents = "wt-12w"), assay = "RNA")
sce2 <- as.SingleCellExperiment(subset(obj, idents = "dare-12w"), assay = "RNA")

metadata(sce1) <- data.frame(colData(sce1))
metadata(sce2) <- data.frame(colData(sce2))

cp1 <- contourPlot(sce1, 
                   dimColnames = c('UMAP_1','UMAP_2'),
                   reducedDim = 'UMAP', lowcol = "white", highcol = "red",
                   bins = 350,
                   subtitle = '',
                   legendLabSize = 18,
                   axisLabSize = 22,
                   titleLabSize = 22,
                   subtitleLabSize = 18,
                   captionLabSize = 18) + 
                                          ggtitle("") + 
                                          theme_classic() + 
                                          theme(panel.grid = element_blank(), 
                                          text = element_blank(), 
                                          axis.ticks = element_blank(), 
                                          legend.title = element_text())

cp2 <- contourPlot(sce2, 
                   dimColnames = c('UMAP_1','UMAP_2'),
                   reducedDim = 'UMAP', lowcol = "white", highcol = "red",
                   bins = 350,
                   subtitle = '',
                   legendLabSize = 18,
                   axisLabSize = 22,
                   titleLabSize = 22,
                   subtitleLabSize = 18,
                   captionLabSize = 18) + 
                                          ggtitle("") + 
                                          theme_classic() + 
                                          theme(panel.grid = element_blank(), 
                                                axis.text.x = element_blank(),
                                                axis.text.y = element_blank(), 
                                                legend.text = element_blank(), 
                                                plot.subtitle = element_blank(),
                                                axis.ticks = element_blank(), 
                                                legend.title = element_text())
cp2$labels$caption <- NULL
cowplot::plot_grid(cp1, cp2,
                   labels = c('Wt','DARE'),
                   ncol = 1, align = "v", label_size = 24)
ggsave("Fig1d.pdf", device = "pdf", width = 25, height = 40, units = "cm", dpi = 600)
