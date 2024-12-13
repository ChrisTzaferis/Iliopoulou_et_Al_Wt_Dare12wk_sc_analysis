library(decoupleR)
library(progeny)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

#----PROGENy model
netMouse <- read.delim("utilityFunctions/progeny_top100_mouse.txt")

obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")

table(obj$sample) #keep dare-12w
Idents(obj) <- obj$sample
obj <- subset(obj, idents = "dare-12w")
Idents(obj) <- obj$annotation
table(obj$annotation)

#----PROGENy model
net <- netMouse

#----Activity inference with Weighted Mean
# Extract the normalized log-transformed counts
obj <- NormalizeData(object = obj, assay = "RNA")
mat <- as.matrix(obj@assays$RNA@data)

# Run wmean with decoupleR
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
tail(acts)
table(acts$source)

#----Visualization
# Extract norm_wmean and store it in pathwayswmean in data
obj[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = obj) <- "pathwayswmean"

# Scale the data
obj <- ScaleData(obj)
obj@assays$pathwayswmean@data <- obj@assays$pathwayswmean@scale.data


#---Plot violin plot---
obj$annotation <- factor(obj$annotation, levels = c("FRC1", "FRC2", 
                                                "Pdgfra_low", "Trophocytes", 
                                                "Telocytes", "Mesothelial", 
                                                "SMCs", "FDCs_MRCs",
                                                "Glial", "Vein_ECs", 
                                                "Pericytes", "LECS", 
                                                "Capillary_ECs", "Cajal"))
Idents(obj) <- obj$annotation

#---superclusters---
new_names <- c("Fibroblasts", "Fibroblasts", 
               "Fibroblasts", "Fibroblasts", 
               "Fibroblasts", "Mesothelial", 
               "SMCs/MFs", "Fibroblasts",
               "Glial", "Endothelial cells",
               "Pericytes", "Endothelial cells", 
               "Endothelial cells", "Cajal") 
names(new_names) <- levels(obj)
obj <- RenameIdents(object = obj, new_names)
obj$supercluster <- Idents(obj)

#---Plot and save---
p1 <- VlnPlot(obj, features = "TNFa", pt.size = 0, group.by = "annotation", 
              cols = c("red3", "red3", 
                       "red3", "red3", 
                       "red3", "grey70", 
                       "orange1", "red3",
                       "forestgreen", "royalblue", 
                       "purple", "royalblue", 
                       "royalblue", "white")) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60)) +
  labs(y="Pathway Activity (z-score)", x="Cluster", title = expression(italic(Tnf^{DARE}) ~ "TNF pathway activity") ) + 
  geom_hline(yintercept = 0, linetype="dashed", color = "black", linewidth=1.2) + NoLegend()

obj$supercluster <- factor(obj$supercluster, levels = c("Fibroblasts", "Mesothelial", "SMCs/MFs", "Glial", 
                                                        "Endothelial cells", "Pericytes", "Cajal"))
p2 <- VlnPlot(obj, features = "TNFa", pt.size = 0, group.by = "supercluster", 
              cols = c("red3", "grey70", "orange1", "forestgreen", "royalblue", "purple", "white"))

#combine plot and legend
final_plot <- plot_grid(
  p1, get_legend(p2), 
  rel_widths = c(2, 0.5),
  ncol = 2
)
final_plot
ggsave("Fig3g.pdf", width = 27, height = 15, dpi = 600, units = "cm")
