library(CellChat)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

pathOI <- "CCL"
pairLR <- extractEnrichedLR(object.list[[1]], signaling = pathOI, geneLR.return = F)
ligandOI <- "Ccl2"
LR.show <- grep(pattern = paste0(ligandOI, "_"), x = pairLR$interaction_name, ignore.case = T, value = T)
sourcesCl <- c("Trophocytes", "Pdgfra_low", "Telocytes")

pdf("Fig4e.pdf", width = 12, height =12)
netVisual_individual(object.list[[1]], signaling = "CCL", pairLR.use = "CCL2_CCR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )

netVisual_individual(object.list[[2]], signaling = "CCL", pairLR.use = "CCL2_CCR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes"), )

netVisual_individual(object.list[[1]], signaling = "CXCL", pairLR.use = "CXCL1_CXCR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes") )

netVisual_individual(object.list[[2]], signaling = "CXCL", pairLR.use = "CXCL1_CXCR2", layout = "circle", 
                     remove.isolate = T, sources.use = c("Trophocytes", "Pdgfra_low", "Telocytes"), )
dev.off()
