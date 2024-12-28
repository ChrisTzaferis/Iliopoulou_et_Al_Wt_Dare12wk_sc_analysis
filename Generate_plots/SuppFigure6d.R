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

obj <- readRDS("../Step9_CellToCellCommunicationAnalysis/Stroma_Myeloid_Lymphoid_obj.RDS")

table(obj$sample) #keep dare-12w
Idents(obj) <- obj$sample
obj <- subset(obj, idents = "dare-12w")
Idents(obj) <- obj$annotation
table(obj$annotation)
obj <- subset(obj, idents = c("FRC1", "FRC2", 
                              "Pdgfra_low", "Trophocytes", 
                              "Telocytes", "Mesothelial", 
                              "SMCs", "FDCs_MRCs",
                              "Glial", "Vein_ECs", 
                              "Pericytes", "LECS", 
                              "Capillary_ECs", "Cajal"), invert=T)

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

#Violin plot
obj$annotation <- factor(obj$annotation, 
                                    levels = c("B_cells","Cycling_B_cells","Plasma_cells", "Naïve_CD4_T_cells", "Naïve_CD8_T_cells",
                                               "Cd8_T_cells", "Memory_T_cells", "Th17_cells","Tregs", "Tcr_gammadelta", "Cycling_T_cells",
                                               "NKT_cells","ILC2s","ILC3s", "Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                                               "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs"
                                    ))

VlnPlot(obj, features = "TNFa", group.by = "annotation", pt.size = 0, 
        cols = c("#CEB2DF", "#E4A149", "#7B78D9", "#B940E1", "#D9A99C", "#D5E2D6", 
                 "#70B3D6", "#D76DC6", "#DC637B", "#70E16E", "#7BE9CE", "#D0E24B",
                 "#D2DE95", "#7C8E66", 
                 "chocolate2","cornflowerblue",
                 "aquamarine4", "darkgoldenrod", "red3", 
                 "darkgreen", "hotpink3","lightpink3", "ivory4")) + 
  theme(legend.position = "none") +
  labs(y="Pathway Activity (z-score)", x="", title = "TNF pathway activity")
ggsave("SFig6d.pdf", device = "pdf", width = 35, height = 15, units = "cm", dpi = 600)
