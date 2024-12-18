library(CellChat)
library(patchwork)

cellchat.NL <- readRDS("cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

levels(cellchat@idents$DARE12)
# 1] "Activated_Mfs"     "B_cells"           "Cajal"             "Capillary_ECs"     "Cd8_T_cells"       "Cdc1"              "Cdc2_a"            "Cdc2_b"           
# [9] "Cycling_B_cells"   "Cycling_T_cells"   "FDCs_MRCs"         "FRC1"              "FRC2"              "Granulocytes"      "ILC2s"             "ILC3s"            
# [17] "Intermediate_Mfs"  "LECS"              "Memory_T_cells"    "Mesothelial"       "Monocytes"         "Naïve_CD4_T_cells" "Naïve_CD8_T_cells" "NKT_cells"        
# [25] "pDCs"              "Pdgfra_low"        "Pericytes"         "Plasma_cells"      "Resident_Mfs"      "SMCs"              "Tcr_gammadelta"    "Telocytes"        
# [33] "Th17_cells"        "Tregs"             "Trophocytes"       "Vein_ECs"

levels(cellchat@idents$WT12)
# [1] "Activated_Mfs"     "B_cells"           "Cajal"             "Capillary_ECs"     "Cd8_T_cells"       "Cdc1"              "Cdc2_a"            "Cdc2_b"           
# [9] "Cycling_B_cells"   "Cycling_T_cells"   "FDCs_MRCs"         "FRC1"              "FRC2"              "Granulocytes"      "ILC2s"             "ILC3s"            
# [17] "Intermediate_Mfs"  "LECS"              "Memory_T_cells"    "Mesothelial"       "Monocytes"         "Naïve_CD4_T_cells" "Naïve_CD8_T_cells" "NKT_cells"        
# [25] "pDCs"              "Pdgfra_low"        "Pericytes"         "Plasma_cells"      "Resident_Mfs"      "SMCs"              "Tcr_gammadelta"    "Telocytes"        
# [33] "Th17_cells"        "Tregs"             "Trophocytes"       "Vein_ECs"     

levels(cellchat@meta$datasets)
# [1] "WT12"   "DARE12"

#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

saveRDS(cellchat, "cellchatCombinedDARE_WT.RDS")