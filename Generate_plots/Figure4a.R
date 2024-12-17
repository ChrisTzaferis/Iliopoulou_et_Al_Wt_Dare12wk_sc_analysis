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

#Differential number of interactions or interaction strength among different cell types
group.cellType <- c("MFs", "B_Plasma cells", "Other", "ECs","Tcells", "DCs", "DCs", "DCs", "B_Plasma cells", 
                    "Tcells", "FB", "FB","FB", "Granulocytes",  "ILCs","ILCs","MFs", "ECs", 
                    "Tcells", "Other", "MFs","Tcells", "Tcells", "Tcells", "DCs", "FB", "Other", 
                    "B_Plasma cells", "MFs", "Other", "Tcells","FB","Tcells", "Tcells", "FB","ECs")

group.cellType <- factor(group.cellType, levels = unique(
  c("MFs", "B_Plasma cells", "Other", "ECs","Tcells", "DCs", "DCs", "DCs", "B_Plasma cells", 
    "Tcells", "FB", "FB","FB", "Granulocytes",  "ILCs","ILCs","MFs", "ECs", 
    "Tcells", "Other", "MFs","Tcells", "Tcells", "Tcells", "DCs", "FB", "Other", 
    "B_Plasma cells", "MFs", "Other", "Tcells","FB","Tcells", "Tcells", "FB","ECs")
))

#Circos plot
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12,  
#                    title.name = paste0("Number of interactions - ", names(object.list)[i]))
# }

pdf("Fig4a.pdf", width = 12, height =12)
netVisual_circle(object.list[[1]]@net$count.merged, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12,  
                 title.name = expression(bolditalic('Tnf'^'+/+')) )

netVisual_circle(object.list[[2]]@net$count.merged, weight.scale = T, label.edge= F, edge.weight.max = weight.max[3], edge.width.max = 12,  
                 title.name = expression(bolditalic('Tnf'^'DARE')) )
dev.off()
