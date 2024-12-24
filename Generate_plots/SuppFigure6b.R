library(Seurat)
library(scCustomize)
library(dplyr)
#---Calculation of DEGs---
#load object
obj.combined <- readRDS("../Step10_Comparison_with_acuteDSS/wt12_HC_dare12_DSSacute_Stroma_integration.RDS")

#cluster_sample column 
obj.combined$cluster_sample <- paste0(obj.combined$annotation, "_", obj.combined$sample)
table(obj.combined$cluster_sample)

ctrH2O <- sort(grep(x = unique(obj.combined$cluster_sample), pattern = "HC", value = T))
expDSS <- sort(grep(x = unique(obj.combined$cluster_sample), pattern = "DSS", value = T))

ctrWT <- sort(grep(x = unique(obj.combined$cluster_sample), pattern = "wt", value = T))
expDARE <- sort(grep(x = unique(obj.combined$cluster_sample), pattern = "dare", value = T))

#calculate DEGs
Idents(obj.combined) <- obj.combined$cluster_sample

#DSS vs H2O
for(i in 1:length(ctrH2O))
{
  degs <- FindMarkers(obj.combined, ident.1 = expDSS[i], ident.2 = ctrH2O[i], assay = "RNA", logfc.threshold = 0.25, min.pct = 0.1, base = exp(1))
  degs$gene <- rownames(degs)
  degs$comparison <- paste0(expDSS[i], "_VS_", ctrH2O[i])
  degs <- degs[which(degs$p_val < 0.01), ]
  degs$cluster <- gsub(x = expDSS[i], pattern = "_DSS_acute", replacement = "")
  write.table(degs, file = paste0("DEGs_", expDSS[i], "_VS_", ctrH2O[i], ".txt"), quote = F, sep = "\t", row.names = F)
}

#DARE vs WT
for(i in 1:length(ctrWT))
{
  degs <- FindMarkers(obj.combined, ident.1 = expDARE[i], ident.2 = ctrWT[i], assay = "RNA", logfc.threshold = 0.25, min.pct = 0.1, base = exp(1))
  degs$gene <- rownames(degs)
  degs$comparison <- paste0(expDARE[i], "_VS_", ctrWT[i])
  degs <- degs[which(degs$p_val < 0.01), ]
  degs$cluster <- gsub(x = expDARE[i], pattern = "_dare-12w", replacement = "")
  write.table(degs, file = paste0("DEGs_", expDARE[i], "_VS_", ctrWT[i], ".txt"), quote = F, sep = "\t", row.names = F)
}

#---Plot venn diagramms---
library(BioVenn)
library(tidyverse)

S1_degs <- read.delim("DEGs_S1_DSS_acute_VS_S1_HC.txt")
S2_degs <- read.delim("DEGs_S2_DSS_acute_VS_S2_HC.txt")
S3_degs <- read.delim("DEGs_S3_DSS_acute_VS_S3_HC.txt")

pdlow <- read.delim("DEGs_Pdgfra_low_dare-12w_VS_Pdgfra_low_wt-12w.txt")
telo <- read.delim("DEGs_Telocytes_dare-12w_VS_Telocytes_wt-12w.txt")
tropho <- read.delim("DEGs_Trophocytes_dare-12w_VS_Trophocytes_wt-12w.txt")


pdf(file = "SFig6b.pdf", width = 9, height = 9)
biovenn <- draw.venn(list_x=NULL, 
                     list_y= pdlow |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S1_degs |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S1 UP", 
                     ytitle = "PDGFRA low UP",
                     z_c = "firebrick", 
                     y_c = "cornflowerblue"
)

biovenn <- draw.venn(list_x=NULL, 
                     list_y= telo |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S2_degs |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S2 UP", 
                     ytitle = "Telocytes UP",
                     z_c = "orange3", 
                     y_c = "cyan3"
)

biovenn <- draw.venn(list_x=NULL, 
                     list_y= tropho |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S3_degs |> filter(avg_logFC > 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S3 UP", 
                     ytitle = "Trophocytes UP",
                     z_c = "yellow3", 
                     y_c = "purple1"
)

#Plotting the downregulated
biovenn <- draw.venn(list_x=NULL, 
                     list_y= pdlow |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S1_degs |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S1 DOWN", 
                     ytitle = "PDGFRA low DOWN",
                     z_c = "firebrick", 
                     y_c = "cornflowerblue"
)

biovenn <- draw.venn(list_x=NULL, 
                     list_y= telo |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S2_degs |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S2 DOWN", 
                     ytitle = "Telocytes DOWN",
                     z_c = "orange3", 
                     y_c = "cyan3"
)

biovenn <- draw.venn(list_x=NULL, 
                     list_y= tropho |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(), 
                     list_z= S3_degs |> filter(avg_logFC < 0.25) |> select(gene) |> as.vector() |> unlist(),
                     subtitle=NULL, 
                     title = NULL, 
                     ztitle = "S3 DOWN", 
                     ytitle = "Tropho DOWN",
                     z_c = "yellow3", 
                     y_c = "purple1"
)
dev.off()