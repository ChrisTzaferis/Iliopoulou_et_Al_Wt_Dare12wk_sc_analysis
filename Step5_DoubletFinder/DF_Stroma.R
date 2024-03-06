library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/wt12_Stroma_PCs_20_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:20, sct = FALSE)
sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
bcmvn_seu <- find.pK(sweep.stats_seu)

bcmvn_seu %>%
  ggplot( aes(x=pK, y=BCmetric, group=1)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=2) +
  theme_bw() +
  ggtitle("pK selection")

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu@meta.data$RNA_snn_res.0.5)
nExp_poi <- round(0.046*nrow(seu@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:20, pN = 0.25, pK = 0.3, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:20, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_453", sct = FALSE)

seu$DF.classifications_0.25_0.3_453 <- gsub(x = seu$DF.classifications_0.25_0.3_453, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.3_453 <- gsub(x = seu$DF.classifications_0.25_0.3_453, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.3_397 <- gsub(x = seu$DF.classifications_0.25_0.3_397, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.3_397 <- gsub(x = seu$DF.classifications_0.25_0.3_397, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.3_453) | as.logical(seu$DF.classifications_0.25_0.3_397)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "wt12_Stroma_PCs_20_0.05_DoubletFinder.RDS")

################################################################################
################################################################################
rm(list=ls(all=TRUE))#remove previous vars, objects

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/dare12_Stroma_PCs_22_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:22, sct = FALSE)
sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
bcmvn_seu <- find.pK(sweep.stats_seu)

bcmvn_seu %>%
  ggplot( aes(x=pK, y=BCmetric, group=1)) +
  geom_line( color="grey") +
  geom_point(shape=21, color="black", fill="#69b3a2", size=2) +
  theme_bw() +
  ggtitle("pK selection")

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seu@meta.data$RNA_snn_res.0.5)           
nExp_poi <- round(0.046*nrow(seu@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:22, pN = 0.25, pK = 0.25, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:22, pN = 0.25, pK = 0.25, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.25_220", sct = FALSE)

seu$DF.classifications_0.25_0.25_220 <- gsub(x = seu$DF.classifications_0.25_0.25_220, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.25_220 <- gsub(x = seu$DF.classifications_0.25_0.25_220, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.25_191 <- gsub(x = seu$DF.classifications_0.25_0.25_191, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.25_191 <- gsub(x = seu$DF.classifications_0.25_0.25_191, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.25_220) | as.logical(seu$DF.classifications_0.25_0.25_191)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "dare12_Stroma_PCs_22_0.05_DoubletFinder.RDS")
