library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)

#Load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
DefaultAssay(obj) <- "RNA"

plot_FP <- function(seu_obj, genes)
{
  g1 <- FeaturePlot(object = seu_obj, features = genes[1], cols = c("grey70", "red"), pt.size = 0.5, order = T) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank()
    )
  g1
  g2 <- FeaturePlot(object = seu_obj, features = genes[2], cols = c("grey70", "red"), pt.size = 0.5, order = T)+ 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank()
    )
  g2
  
  combined_plot <- (g1 | g2) + plot_annotation(theme = theme(plot.background = element_rect(color = "black", size = 1.5)))
  return(combined_plot)
}

pdf(file = "SFig5a.pdf", width = 10, height = 5)
plot_FP(seu_obj = obj, genes = c("Pecam1", "Cdh5"))
plot_FP(seu_obj = obj, genes = c("Pdgfra", "Pdpn"))
plot_FP(seu_obj = obj, genes = c("Vtn", "Pdgfrb"))
plot_FP(seu_obj = obj, genes = c("Wt1", "Upk3b"))
plot_FP(seu_obj = obj, genes = c("Kit", "Ano1"))
plot_FP(seu_obj = obj, genes = c("Cnn1", "Myh11"))
plot_FP(seu_obj = obj, genes = c("Sox10", "S100b"))
dev.off()
