library(Seurat)
library(dplyr)
library(tidyr)
library(UCell)
library(scCustomize)
library(ggplot2)
library(nichenetr)
library(dittoSeq)
library(cowplot)

obj1 <- readRDS("../Step8_Preparation_of_human_dataset/myeloid_subclusters_annotated.RDS") 
obj2 <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")

#Markers genes for all clusters
hMarkers <- FindAllMarkers(obj1, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.25, only.pos = T, return.thresh = 0.01, base = exp(1)) 
#-------------------------------------------------------------------------------

#Human to mouse conversion (using nichnetR function) and selection of top50 markers per cluster
top50 <- hMarkers |> 
  group_by(cluster) |> 
  slice_max(avg_logFC, n = 50)

unique(top50$gene) #332 from 350 unique

top50 <- top50[, c("gene", "cluster")]
top50M <- top50 |> 
  mutate(geneMouse = convert_human_to_mouse_symbols(gene)) |>
  relocate(gene, geneMouse) |>
  dplyr::select(geneMouse, cluster)
top50M <- na.omit(top50M)
top50M |> group_by(cluster) |> summarise(n=n())

#Create a list for UCell
ll <- split(top50M, f = top50M$cluster)
ll <- lapply(X = ll, FUN = function(x){
  x$cluster <- NULL
  x <- as.vector(x$geneMouse)
})

#calculate signatures
obj2 <- AddModuleScore_UCell(obj2, features = ll, assay = "RNA", slot = "data")
signature.names <- paste0(names(ll), "_UCell")

obj2$cluster_sample <- paste0(obj2$annotation, "_", obj2$sample)

#Plots for the final figures
cellsOfInterest <- WhichCells(object = obj2, idents = c("Monocytes", "Intermediate_Mfs", "Resident_Mfs", "Activated_Mfs"))
p1 <- dittoPlot(obj2, "Resident_Macs_UCell", group.by = "cluster_sample", cells.use = cellsOfInterest, 
          color.panel = rep(c("red", "grey70"), 4), 
          plots = c("vlnplot", "boxplot", "jitter"),
          #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
          
          # change the color and size of jitter points
          jitter.color = "black", jitter.size = 0.6,
          
          # change the outline color and width, and remove the fill of boxplots
          boxplot.color = "royalblue", boxplot.width = 0.5,
          boxplot.fill = FALSE, boxplot.show.outliers = T,
          
          # violin plot widths
          vlnplot.scaling = "width", vlnplot.lineweight = 0.8, 
          
          #labels
          x.reorder = c(5, 6, 3, 4, 7, 8, 1, 2), theme = theme_bw(), legend.show = F,
          main = "Human resident MFs signature (HC/CD)", xlab = "", ylab = "UCell signature score"
) + theme(text = element_text(size = 16))

p2 <- dittoPlot(obj2, "Inflammatory_Macs_UCell", group.by = "cluster_sample", cells.use = cellsOfInterest, 
          color.panel = rep(c("red", "grey70"), 4), 
          plots = c("vlnplot", "boxplot", "jitter"),
          #plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
          
          # change the color and size of jitter points
          jitter.color = "black", jitter.size = 0.6,
          
          # change the outline color and width, and remove the fill of boxplots
          boxplot.color = "royalblue", boxplot.width = 0.5,
          boxplot.fill = FALSE, boxplot.show.outliers = T,
          
          # violin plot widths
          vlnplot.scaling = "width", vlnplot.lineweight = 0.8, 
          
          #labels
          x.reorder = c(5, 6, 3, 4, 7, 8, 1, 2), theme = theme_bw(), legend.show = F,
          main = "Human pro-inflammatory MFs signature (HC/CD)", xlab = "", ylab = "UCell signature score"
) + theme(text = element_text(size = 16), axis.text.x = element_blank()) 

#get legend from a violin plot
p3 <- VlnPlot( obj2, features = "Inflammatory_Macs_UCell", cols = c("red", "grey70"), group.by = "sample") + 
  scale_fill_manual(
    values = c(
      "wt-12w" = "grey70",
      "dare-12w" = "red"
    ),
    labels = c(
      "wt-12w" = expression(italic('Tnf'^'+/+')),
      "dare-12w" = expression(italic('Tnf'^'DARE'))
    )
  ) + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    legend.text.align = 0,
    text = element_text(size = 22)
  )

# Combine plots with the legend
aligned_plots <- plot_grid(
  p2, p1,
  align = "v",
  axis = "l",   
  ncol = 1, rel_heights = c(1, 1.4)
)

# Add legend next to plots
final_plot <- plot_grid(
  aligned_plots, get_legend(p3), 
  rel_widths = c(2, 0.5),
  ncol = 2
)
final_plot
ggsave("Figure2k.pdf", width = 29, height = 29, dpi = 600, units = "cm")

