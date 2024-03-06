library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)

# Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/wt12_Lymphoid_PCs_15_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:15, sct = FALSE)
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
seu <- doubletFinder_v3(seu, PCs = 1:15, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:15, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.29_539", sct = FALSE)

seu$DF.classifications_0.25_0.29_539 <- gsub(x = seu$DF.classifications_0.25_0.29_539, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.29_539 <- gsub(x = seu$DF.classifications_0.25_0.29_539, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.29_471 <- gsub(x = seu$DF.classifications_0.25_0.29_471, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.29_471 <- gsub(x = seu$DF.classifications_0.25_0.29_471, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.29_539) | as.logical(seu$DF.classifications_0.25_0.29_471)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "wt12_Lymphoid_PCs_15_0.05_DoubletFinder.RDS")

################################################################################
################################################################################
rm(list=ls(all=TRUE))#remove previous vars, objects

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- readRDS("../Step4_Analysis_of_Stroma_Lymphoid_Myeloid/dare12_Lymphoid_PCs_23_0.5.RDS")

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_seu <- paramSweep_v3(seu, PCs = 1:23, sct = FALSE)
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
seu <- doubletFinder_v3(seu, PCs = 1:23, pN = 0.25, pK = 0.05, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:23, pN = 0.25, pK = 0.05, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.05_393", sct = FALSE)

seu$DF.classifications_0.25_0.05_393 <- gsub(x = seu$DF.classifications_0.25_0.05_393, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.05_393 <- gsub(x = seu$DF.classifications_0.25_0.05_393, pattern = "Singlet", replacement = FALSE)
seu$DF.classifications_0.25_0.05_316 <- gsub(x = seu$DF.classifications_0.25_0.05_316, pattern = "Doublet", replacement = TRUE)
seu$DF.classifications_0.25_0.05_316 <- gsub(x = seu$DF.classifications_0.25_0.05_316, pattern = "Singlet", replacement = FALSE)

seu$DoubletFinderClass <- as.logical(seu$DF.classifications_0.25_0.05_316) | as.logical(seu$DF.classifications_0.25_0.05_316)

DimPlot(seu, group.by = "DoubletFinderClass", split.by = "DoubletFinderClass")
saveRDS(seu, "dare12_Lymphoid_PCs_23_0.05_DoubletFinder.RDS")
