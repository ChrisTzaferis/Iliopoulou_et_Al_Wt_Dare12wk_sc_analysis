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

#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = c("grey70", "red2")) 
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("grey70", "red2"))
gg1 + gg2 & theme(text = element_text(size = 12)) & scale_x_discrete(labels=c("WT12" = expression(italic("Tnf"^"+/+")), "DARE12" = expression(italic("Tnf"^"DARE")) ))
ggsave("SFig7a.pdf", device = "pdf", width = 12, height = 10, units = "cm", dpi = 600)
