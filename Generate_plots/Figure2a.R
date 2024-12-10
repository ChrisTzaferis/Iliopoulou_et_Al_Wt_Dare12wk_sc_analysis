library(Seurat)
library(scCustomize)
library(ggplot2)
library(cowplot)
library(patchwork)

#---load object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                 "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))

obj_list <- SplitObject(obj, split.by = "sample")

#---UMAP plot for WT---
umap1 <- DimPlot_scCustom(obj_list[[1]], group.by = "annotation", colors_use = c("chocolate2","cornflowerblue",
                                                                        "aquamarine4", "darkgoldenrod", "red3", 
                                                                        "darkgreen", "hotpink3","lightpink3", "ivory4"), 
                 figure_plot = T, 
                 split_seurat = F,
                 label = T, 
                 pt.size = 0.8, 
                 repel = T
) & theme(legend.position = "none", 
          line = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25))

#---UMAP for DARE---
umap2 <- DimPlot_scCustom(obj_list[[2]], group.by = "annotation", colors_use = c("chocolate2","cornflowerblue",
                                                                                 "aquamarine4", "darkgoldenrod", "red3", 
                                                                                 "darkgreen", "hotpink3","lightpink3", "ivory4"), 
                          figure_plot = T, 
                          split_seurat = F,
                          label = T, 
                          pt.size = 0.8, 
                          repel = T
) & theme(legend.position = "none", 
          line = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 25))

#UMAP plot to get the color legend
umap3 <- DimPlot_scCustom(obj, group.by = "annotation", colors_use = c("chocolate2","cornflowerblue",
                                                                       "aquamarine4", "darkgoldenrod", "red3", 
                                                                       "darkgreen", "hotpink3","lightpink3", "ivory4")
)
col_legend <- get_legend(umap3)

#---customise WT umap--- 
plot1 <- umap1[[2]] 
plot2 <- umap1[[1]] + ggtitle(expression(italic('Tnf'^'+/+')))

# Create an empty placeholder for the top part of the umap axes
empty_top <- plot_grid(NULL, ncol = 1)

# Create the axes with an empty top area
small_plot_with_space <- plot_grid(
  empty_top, plot1, 
  nrow = 2, 
  rel_heights = c(3, 1) # Empty space (3/4) above the umap axes (1/4)
)

# Combine the umap axes (on the left) with the umap scatter (on the right)
g1 <- plot_grid(
  small_plot_with_space, plot2, 
  ncol = 2, 
  rel_widths = c(1, 3) # Small plot takes 1/4 of the width
)

#---customise DARE umap--- 
plot3 <- umap2[[2]] 
plot4 <- umap2[[1]] + ggtitle(expression(italic('Tnf'^'DARE')))

# Create an empty placeholder for the top part of the umap axes
empty_top2 <- plot_grid(NULL, ncol = 1)

# Create the axes with an empty top area
small_plot_with_space2 <- plot_grid(
  empty_top2, plot3, 
  nrow = 2, 
  rel_heights = c(3, 1) # Empty space (3/4) above the umap axes (1/4)
)

# Combine the umap axes (on the left) with the umap scatter (on the right)
g2 <- plot_grid(
  small_plot_with_space2, plot4, 
  ncol = 2, 
  rel_widths = c(1, 3) # Small plot takes 1/4 of the width
)

#Plot with legend and save
plot_grid(
  g1+g2,
  col_legend,                         
  ncol = 2,                            
  rel_widths = c(5, 1)
)
ggsave("Fig2a.pdf", device = "pdf", width = 39, height = 17, units = "cm", dpi = 600)
