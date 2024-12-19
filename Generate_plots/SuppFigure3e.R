library(Seurat)
library(ggplot2)
library(slingshot)
library(dittoSeq)
library(dplyr)
library(viridis)
library(Hmisc)

#load object and keep only clusters of macrophages
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
DimPlot(obj, label = T)

obj <- subset(obj, idents=c("Monocytes", "Intermediate_Mfs", "Resident_Mfs", "Activated_Mfs"))
obj$sample <- factor(obj$sample, levels = c("wt-12w", "dare-12w"))
obj$annotation <- factor(obj$annotation)

#downsample obj
obj_list <- SplitObject(obj, split.by = "sample")
cell_ids_wt <- WhichCells(obj_list[[1]])
dare_downsampled <- subset(obj_list[[2]], cells = sample(Cells(obj_list[[2]]), 1222)) #length(cell_ids_wt)=1222
obj <- subset(obj, cells = c(WhichCells(obj_list[[1]]), WhichCells(dare_downsampled)))

#Slingshot input
dimred1 <- obj@reductions$pca@cell.embeddings[, 1:20]
dimred2 <- obj@reductions$umap@cell.embeddings
clustering <- obj$annotation

#Run default Slingshot
set.seed(1)
dimred <- dimred2
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        omega = T,
                        #end.clus = c(), #define the end of the trajectory
                        start.clus = "Monocytes") #define where to start the trajectories

#Exploring the results
lineages

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0, allow.breaks = FALSE, shrink = 0.99)
curves

dittoDimPlot(obj, var = "sample", add.trajectory.lineages = slingLineages(lineages), trajectory.cluster.meta = "annotation", 
             do.label = F, split.by = "sample", show.others = F, size = 1.2)
slingLineages(lineages)

pal <- dittoColors()[1:4]
plot(dimred, asp = 1, pch = 16, col = pal[clustering])
lines(SlingshotDataSet(curves), lwd = 3, col = 'black')

#---Pseudotime values---
pto <- as.data.frame(curves@assays@data$pseudotime[, c(1, 2)])
pto$Cell_id <- rownames(pto)
cluster_info <- as.data.frame(obj@meta.data[, c("annotation", "sample")])
cluster_info$Cell_id <- rownames(cluster_info)

#dataframe with pseudotime, clusters, sample info for the two lineages
cluster_info <- left_join(cluster_info, pto)

#Check NAs pseudotime values in each lineage 
nasInLin1 <- cluster_info[which(is.na(cluster_info$Lineage1)), ] 
table(nasInLin1$annotation)

nasInLin2 <- cluster_info[which(is.na(cluster_info$Lineage2)), ] 
table(nasInLin2$annotation)

#---Remove NAs in lineage1 and divide time in 5 bins---
notNasInLin1 <- cluster_info[which(is.na(cluster_info$Lineage1) == FALSE), ] 
notNasInLin1$bin <- as.numeric(cut2(notNasInLin1$Lineage1, g=5))
notNasInLin1$bin <- paste0("Bin_", notNasInLin1$bin)

#calculate values for the first barplot
final_df1 <- notNasInLin1
final_df1 <- final_df1 %>% group_by(bin, sample) %>%
  summarise(N_cells=n())

total_cells_in_df1 <- as.data.frame(table(notNasInLin1$bin))
colnames(total_cells_in_df1) <- c("bin", "Total_cells_in_bin")

final_df1 <- left_join(final_df1, total_cells_in_df1)
final_df1$Percentage <- (final_df1$N_cells/final_df1$Total_cells_in_bin)*100
final_df1$sample <- factor(final_df1$sample, levels = c("wt-12w", "dare-12w"))
final_df1$bin <- factor(final_df1$bin, levels = c("Bin_1", "Bin_2", "Bin_3", "Bin_4", "Bin_5"))

#first barplot
g1 <- ggplot(final_df1, aes(y=Percentage, fill=sample, x=sample )) +
  theme_bw() + 
  facet_wrap(~bin, ncol = 5) + 
  geom_bar(stat = "identity" ) +
  scale_fill_manual(values = c("darkgreen", "red4"), 
                    labels = c("wt-12w"= expression('Tnf'^'+/+'), "dare-12w" = expression('Tnf'^'DARE'))) +
  theme(legend.position="right", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0,
        legend.title = element_blank(),
        text = element_text(size = 25)) +
  labs(x="Sample", y="Percentage of cells", title = "Monocyte to Activated Macrophage lineage")


#---Remove NAs in lineage2---
notNasInLin2 <- cluster_info[which(is.na(cluster_info$Lineage2) == FALSE), ] 
notNasInLin2$bin <- as.numeric(cut2(notNasInLin2$Lineage2, g=5))
notNasInLin2$bin <- paste0("Bin_", notNasInLin2$bin)

final_df2 <- notNasInLin2
final_df2 <- final_df2 %>% group_by(bin, sample) %>%
  summarise("N_cells" = n())
total_cells_in_df2 <- as.data.frame(table(notNasInLin2$bin))
colnames(total_cells_in_df2) <- c("bin", "Total_cells_in_bin")
final_df2 <- left_join(final_df2, total_cells_in_df2)
final_df2$Percentage <- (final_df2$N_cells/final_df2$Total_cells_in_bin)*100
final_df2$sample <- factor(final_df2$sample, levels = c("wt-12w", "dare-12w"))
final_df2$bin <- factor(final_df2$bin, levels = c("Bin_1", "Bin_2", "Bin_3", "Bin_4", "Bin_5"))

#second barplot
g2 <- ggplot(final_df2, aes(y=Percentage, fill=sample, x=sample )) +
  theme_bw() + 
  facet_wrap(~bin, ncol = 5) + 
  geom_bar(stat = "identity" ) +
  scale_fill_manual(values = c("darkgreen", "red4"), 
                    labels = c("wt-12w"= expression('Tnf'^'+/+'), "dare-12w" = expression('Tnf'^'DARE'))) +
  theme(legend.position="right", 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0,
        legend.title = element_blank(),
        text = element_text(size = 25)) +
  labs(x="Sample", y="Percentage of cells", title = "Monocyte to Resident Macrophage lineage")

#combine plots and save
g2/g1
ggsave("SFig3e.pdf", width = 29, height = 29, dpi = 600, units = "cm")
