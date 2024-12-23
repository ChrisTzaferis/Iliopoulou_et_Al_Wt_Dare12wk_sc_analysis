library(ggplot2)
library(dplyr)
library(tidyverse)
library(Seurat)
library(decoupleR)

#Load object and keep only fibroblast clusters
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)
Idents(obj) <- obj$cluster_sample
obj <- subset(obj, idents = c("Trophocytes_wt-12w", "Trophocytes_dare-12w",
                              "Pdgfra_low_wt-12w", "Pdgfra_low_dare-12w",
                              "Telocytes_wt-12w", "Telocytes_dare-12w"))
table(obj$cluster_sample)
#----PROGENy model
net <- read.delim("utilityFunctions/progeny_top100_mouse.txt")

#----Activity inference with Weighted Mean
# Extract the normalized log-transformed counts
obj <- NormalizeData(object = obj, assay = "RNA")
mat <- as.matrix(obj@assays$RNA@data)

# Run wmean with decoupleR
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
tail(acts)
table(acts$source)

# Extract norm_wmean and store it in pathwayswmean in data
obj[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Scale the data
DefaultAssay(obj) <- "pathwayswmean"
obj <- ScaleData(obj)
#obj@assays$pathwayswmean@data <- obj@assays$pathwayswmean@scale.data

#Get mean activity values
df <- t(as.matrix(obj@assays$pathwayswmean@scale.data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(obj)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

#Get TNF and TGFbeta
tnf_tgfb_df <- df |> filter(source %in% c("TNFa", "TGFb"))

#convert to wide format and then to matrix
htable <- tnf_tgfb_df %>%
  pivot_wider(names_from = source, values_from = mean)
hmat <- as.matrix(htable[, -1])
rownames(hmat) <- htable$cluster

#plot and save heatmap
palette_length = 20
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-1.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1.5, length.out=floor(palette_length/2)))


pdf("SFig5h.pdf", width = 10, height = 20)
pheatmap::pheatmap(hmat[c(2, 1, 3, 4, 5, 6), c(2,1)], 
                   color = my_color,
                   breaks = my_breaks,
                   scale = "none", 
                   cellwidth = 120, 
                   cellheight = 120,
                   fontsize = 34, 
                   legend_breaks = c(-1.5, 0, 1.5),
                   cluster_rows = F, 
                   cluster_cols = F)
dev.off()
