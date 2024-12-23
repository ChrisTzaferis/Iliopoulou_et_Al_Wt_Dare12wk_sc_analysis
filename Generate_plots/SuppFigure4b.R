source("utilityFunctions/scGSVA_functions.R")

#load myeloid 
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)

#run and save ssGSEA
seu_REAC_myel <- performSCGSVA(seu_obj = obj, geneSetType = "REACTOME", calcMethod = "ssgsea")
#saveRDS(seu_REAC_myel, "myeloid_ssGSEA_REACTOME.RDS")  #(optional)

#plot heatmap with top5
plotTopSD(seu_obj = seu_REAC_myel, groupCluster = "cluster_sample", assayName = "ssGSEA", showAll = F, nToPlot = 5, filenameOut = "SFig4b.pdf")