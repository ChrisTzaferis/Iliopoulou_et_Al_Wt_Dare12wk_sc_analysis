library(decoupleR)
library(progeny)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(Seurat)
library(cowplot)

#----PROGENy model
netHuman <- read.delim("utilityFunctions/progeny_top100_human.txt")

obj <- readRDS("../Step8_Preparation_of_human_dataset/TI_Stroma_Human.RDS")

table(obj$Type) #keep Infl
Idents(obj) <- obj$Type
obj <- subset(obj, idents = "Infl")
Idents(obj) <- obj$Celltype
table(obj$Celltype)

#----PROGENy model
net <- netHuman

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
obj$Celltype <- factor(obj$Celltype, levels = c("Activated fibroblasts CCL19 ADAMADEC1", "Endothelial cells DARC", 
                                                "Fibroblasts ADAMDEC1", "Endothelial cells CD36", 
                                                "Fibroblasts SFRP2 SLPI", "Fibroblasts SMOC2 PTGIS", 
                                                "Fibroblasts NPY SLITRK6", "Endothelial cells CA4 CD36",
                                                "Lymphatics", "Glial cells", 
                                                "Pericytes HIGD1B STEAP4", "Fibroblasts KCNN3 LY6H", 
                                                "Pericytes RERGL NTRK2", "Endothelial cells LTC4S SEMA3G", 
                                                "Myofibroblasts HHIP NPNT", "Myofibroblasts GREM1 GREM2"))
Idents(obj) <- obj$Celltype

#---superclusters---
new_names <- c("Fibroblasts", "Endothelial cells", 
               "Fibroblasts", "Endothelial cells", 
               "Fibroblasts", "Fibroblasts", 
               "Fibroblasts", "Endothelial cells",
               "Endothelial cells", "Glial", 
               "Pericytes", "Fibroblasts", 
               "Pericytes", "Endothelial cells", 
               "SMCs/MFs", "SMCs/MFs") 
names(new_names) <- levels(obj)
obj <- RenameIdents(object = obj, new_names)
obj$supercluster <- Idents(obj)

#---Plot and save---
p1 <- VlnPlot(obj, features = "TNFa", pt.size = 0, group.by = "Celltype", 
              cols = c("red3", "royalblue", 
                       "red3", "royalblue", 
                       "red3", "red3", 
                       "red3", "royalblue",
                       "royalblue", "forestgreen", 
                       "purple", "red3", 
                       "purple", "royalblue", 
                       "orange1", "orange1")) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 60)) +
  labs(y="Pathway Activity (z-score)", x="Cluster", title = "Human CD TNF pathway activity") + 
  geom_hline(yintercept = 0, linetype="dashed", color = "black", linewidth=1.2) + NoLegend()

obj$supercluster <- factor(obj$supercluster, levels = c("Fibroblasts", "SMCs/MFs", "Glial", "Endothelial cells", "Pericytes"))
p2 <- VlnPlot(obj, features = "TNFa", pt.size = 0, group.by = "supercluster", 
              cols = c("red3", "orange1", "forestgreen", "royalblue", "purple"))

#combine plot and legend
final_plot <- plot_grid(
  p1, get_legend(p2), 
  rel_widths = c(2, 0.5),
  ncol = 2
)
final_plot
ggsave("Fig3h.pdf", width = 29, height = 15, dpi = 600, units = "cm")
