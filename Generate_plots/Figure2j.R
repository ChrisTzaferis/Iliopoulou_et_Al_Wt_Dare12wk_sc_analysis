#load function
source("utilityFunctions/heatmapForMetascape.R")

#file produced from metascape
metafile <- read.csv("../Metascape_files/fig2j_GO_AllLists.csv")

#Selected terms
selection1 <- c("positive regulation of cell migration", "Trafficking and processing of endosomal TLR", 
                "endocytosis", "response to wounding", "collagen metabolic process", "inflammatory response", 
                "innate immune response", "chemoattractant activity", "Oxidative stress response", 
                "HIF-1 signaling pathway - Mus musculus (house mouse)", 
                "Protein processing in endoplasmic reticulum - Mus musculus (house mouse)", 
                "chaperone cofactor-dependent protein refolding", "ubiquitin protein ligase binding", 
                "positive regulation of monocyte differentiation")

#selected clusters
selection2 <- c("Activated_Mfs", "Resident_Mfs")

#print heatmap
metascapeToHeatmap(meta_table = metafile, selectedTerms = selection1, selectedClusters = selection2, maxCol = 12, 
                   pdfName = "Fig2j.pdf", pdfWidth = 12, pdfHeight = 12, clusterRows = F, clusterCols = F)
