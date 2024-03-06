library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/wt12_Myeloid_PCs_14_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:14, sct = FALSE)
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
seu <- doubletFinder_v3(seu, PCs = 1:14, pN = 0.25, pK = 0.2, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:14, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_149", sct = FALSE)

seu$DF.classifications_0.25_0.2_149 <- gsub(x = seu$DF.classifications_0.25_0.2_149, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.2_149 <- gsub(x = seu$DF.classifications_0.25_0.2_149, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.2_134 <- gsub(x = seu$DF.classifications_0.25_0.2_134, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.2_134 <- gsub(x = seu$DF.classifications_0.25_0.2_134, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.2_134) | as.logical(seu$DF.classifications_0.25_0.2_149)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "wt12_Myeloid_PCs_14_0.5_DoubletFinder.RDS")

################################################################################
################################################################################
rm(list=ls(all=TRUE))#remove previous vars, objects

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/dare12_Myeloid_PCs_21_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:21, sct = FALSE)
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
seu <- doubletFinder_v3(seu, PCs = 1:21, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:21, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_329", sct = FALSE)

seu$DF.classifications_0.25_0.03_329 <- gsub(x = seu$DF.classifications_0.25_0.03_329, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.03_329 <- gsub(x = seu$DF.classifications_0.25_0.03_329, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.03_284 <- gsub(x = seu$DF.classifications_0.25_0.03_284, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.03_284 <- gsub(x = seu$DF.classifications_0.25_0.03_284, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.03_329) | as.logical(seu$DF.classifications_0.25_0.03_284)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "dare12_Myeloid_PCs_21_0.5_DoubletFinder.RDS")
