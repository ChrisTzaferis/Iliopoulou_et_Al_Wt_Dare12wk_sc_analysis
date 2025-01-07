library(CellChat)
library(tidyverse)
library(ggplot2)
library(patchwork)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

#Plot L-R violin plots for CXCL, CCL
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT12", "DARE12")) # set factor level

a <- plotGeneExpression(cellchat, signaling = c("CCL", "CXCL"), split.by = "datasets", 
                   features = c("Ccl2", "Ccr2", "Ackr1", 
                                "Ackr2", "Cxcl1", "Cxcr2"),
                   colors.ggplot = T, color.use = c("grey70", "red2"))

pdf(file = "SFig11c.pdf", width = 15, height = 20)
a & scale_fill_manual(
  values = c(
    "WT12" = "grey70",
    "DARE12" = "red2"
  ),
  labels = c(
    "WT12" = expression(italic('Tnf'^'+/+')),
    "DARE12" = expression(italic('Tnf'^'DARE'))
  )
) & theme(legend.text.align = 0, text = element_text(size = 20))
dev.off()
