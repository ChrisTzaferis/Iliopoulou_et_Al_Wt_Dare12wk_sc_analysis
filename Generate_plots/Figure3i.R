#load function
source("utilityFunctions/heatmapForMetascape.R")

metafile <- read.csv("../Metascape_files/fig3i_GO_AllLists.csv")
metafile$GeneList <- gsub(pattern = "dare-12w__Fibroblasts_M_vs_wt-12w__Fibroblasts_M_DOWN", replacement = "Mouse_DOWN", x = metafile$GeneList)
metafile$GeneList <- gsub(pattern = "dare-12w__Fibroblasts_M_vs_wt-12w__Fibroblasts_M_UP", replacement = "Mouse_UP", x = metafile$GeneList)
metafile$GeneList <- gsub(pattern = "Infl__Fibroblasts_H_vs_Heal__Fibroblasts_H_DOWN", replacement = "Human_DOWN", x = metafile$GeneList)
metafile$GeneList <- gsub(pattern = "Infl__Fibroblasts_H_vs_Heal__Fibroblasts_H_UP", replacement = "Human_UP", x = metafile$GeneList)

#Selected terms
selection1 <- c("inflammatory response", "cytokine-mediated signaling pathway", "TNF signaling pathway - Mus musculus (house mouse)", 
                "cell activation", "chemokine activity", "response to lipopolysaccharide", "Comprehensive IL-17A signaling", 
                "Elastic fibre formation", "MAPK family signaling cascades", 
                "response to oxidative stress", "regulation of BMP signaling pathway", "smooth muscle cell differentiation", "regulation of sterol transport"
)

#selected clusters
selection2 <- c("Human_UP", "Mouse_UP", "Human_DOWN", "Mouse_DOWN")

#print heatmap
metascapeToHeatmap(meta_table = metafile, selectedTerms = selection1, selectedClusters = selection2, maxCol = 12, 
                   pdfName = "Fig3i.pdf", pdfWidth = 8, pdfHeight = 12, clusterRows = F, clusterCols = F)
