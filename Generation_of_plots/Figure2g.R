library(pheatmap)
library(dplyr)

#file produced from metascape
meta_table <- read.csv("../Metascape_files/fig2g_FINAL_GO.csv")

#selected terms
selection <- c("positive regulation of cell migration", "Trafficking and processing of endosomal TLR", 
               "endocytosis", "response to wounding", "collagen metabolic process", "inflammatory response", 
               "innate immune response", "chemoattractant activity", "Oxidative stress response", 
               "HIF-1 signaling pathway - Mus musculus (house mouse)", 
               "Protein processing in endoplasmic reticulum - Mus musculus (house mouse)", 
               "chaperone cofactor-dependent protein refolding", "ubiquitin protein ligase binding", 
               "positive regulation of monocyte differentiation")

#keep selected terms and -log10Pval columns
meta_table_filtered <- meta_table[which(meta_table$Description %in% selection), ]
heatmap_table <- meta_table_filtered[, c(3, 4, 11)]
rownames(heatmap_table) <- heatmap_table[, 3]
heatmap_table <- as.matrix(heatmap_table[, -3])
heatmap_table <- heatmap_table*(-1)
colnames(heatmap_table) <- c("Activated_Mfs", "Resident_Mfs")

#print heatmap
pdf(file = "Fig2g.pdf", height = 15, width = 8)
pheatmap(heatmap_table[selection, ], cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(c("white", "lightpink","red", "red3", "red4"))(50))
dev.off()
