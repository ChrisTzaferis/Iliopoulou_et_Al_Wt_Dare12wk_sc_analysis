library(CellChat)
library(tidyverse)
library(ggplot2)
library(patchwork)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
levels(cellchat@meta$annotation)
# [1] "Activated_Mfs"     "B_cells"           "Cajal"             "Capillary_ECs"     "Cd8_T_cells"       "Cdc1"              "Cdc2_a"            "Cdc2_b"            "Cycling_B_cells"  
# [10] "Cycling_T_cells"   "FDCs_MRCs"         "FRC1"              "FRC2"              "Granulocytes"      "ILC2s"             "ILC3s"             "Intermediate_Mfs"  "LECS"             
# [19] "Memory_T_cells"    "Mesothelial"       "Monocytes"         "Naïve_CD4_T_cells" "Naïve_CD8_T_cells" "NKT_cells"         "pDCs"              "Pdgfra_low"        "Pericytes"        
# [28] "Plasma_cells"      "Resident_Mfs"      "SMCs"              "Tcr_gammadelta"    "Telocytes"         "Th17_cells"        "Tregs"             "Trophocytes"       "Vein_ECs" 

cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

#The first three pages are for WT, while the rest are for DARE
pdf(file = "SFig11b.pdf", width = 10, height = 10)
netVisual_individual(object.list[[1]], signaling = "CCL", pairLR.use = "CCL2_ACKR1", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") ) 

netVisual_individual(object.list[[1]], signaling = "CCL", pairLR.use = "CCL2_ACKR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )

netVisual_individual(object.list[[1]], signaling = "CXCL", pairLR.use = "CXCL1_ACKR1", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )


netVisual_individual(object.list[[2]], signaling = "CCL", pairLR.use = "CCL2_ACKR1", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") ) 

netVisual_individual(object.list[[2]], signaling = "CCL", pairLR.use = "CCL2_ACKR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )

netVisual_individual(object.list[[2]], signaling = "CXCL", pairLR.use = "CXCL1_ACKR1", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )
dev.off()