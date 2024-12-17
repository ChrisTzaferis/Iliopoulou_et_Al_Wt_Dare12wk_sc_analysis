library(CellChat)
library(ggplot2)
library(patchwork)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

#save scatter plots
p1 <- gg[[1]] + theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 25)) + xlim(0, 100) + ylim(0, 80) + ggtitle(expression(bolditalic('Tnf'^'+/+')))
p2 <- gg[[2]] + theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 25)) + xlim(0, 100) + ylim(0, 80) + ggtitle(expression(bolditalic('Tnf'^'DARE')))

p1+p2
ggsave(filename = "Fig4b.pdf", device = "pdf", width = 32, height = 20, units = "cm", dpi = 600)
