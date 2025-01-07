library(CellChat)
library(ggplot2)

#load objects
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

#Signaling pathways with FDCs/MRCs as source 
pdf("SFig9a.pdf", width = 15, height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = 11, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[1]))

netVisual_chord_gene(object.list[[2]], sources.use = 11, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[2]))
dev.off()

#Signaling pathways with FRC1 as source 
pdf("SFig9b.pdf", width = 15, height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = 12, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[1]))

netVisual_chord_gene(object.list[[2]], sources.use = 12, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[2]))
dev.off()

#Signaling pathways with FRC2 as source 
pdf("SFig9c.pdf", width = 15, height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = 13, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[1]))

netVisual_chord_gene(object.list[[2]], sources.use = 13, targets.use = c(11, 12, 13), slot.name = "netP", small.gap = 3,
                     title.name = paste0("", names(object.list)[2]))
dev.off()