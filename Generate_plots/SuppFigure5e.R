source("utilityFunctions/scGSVA_functions.R")

#load stroma 
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)

#run and save ssGSEA 
seu_HM <- performSCGSVA(seu_obj = obj, geneSetType = "HALLMARK", calcMethod = "ssgsea")
#saveRDS(seu_HM, "stroma_ssGSEA_HALLM.RDS") #(optional)

#plot heatmap with top5
plotTopSD(seu_obj = seu_HM, groupCluster = "cluster_sample", identsToExclude = "Glial_wt-12w", assayName = "ssGSEA", showAll = F, 
          nToPlot = 5, filenameOut = "SFig5e.pdf")
