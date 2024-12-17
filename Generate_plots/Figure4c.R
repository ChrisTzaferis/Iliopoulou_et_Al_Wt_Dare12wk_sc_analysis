library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(Seurat)
library(scales)
library(patchwork)
library(decoupleR)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

#DARE12 weeks
#1)---receiver barplot---
dare_table <- subsetCommunication(cellchat.LS, slot.name = "netP", signaling = "TNF")
dare_sum <- dare_table %>% 
  group_by(target, pathway_name) %>% 
  summarise(prob = sum(prob))
df_label_colors <- data.frame(target=levels(gg[[1]]$data$labels),
                              colors=gg[[1]]$plot_env$color.use)
dare_sum <- left_join(df_label_colors, dare_sum)
dare_sum$prob[is.na(dare_sum$prob)] <- 0
dare_sum$target <- factor(dare_sum$target, levels = rev(
  c("B_cells", "Cycling_B_cells", "Plasma_cells", "ILC2s", "ILC3s", "Granulocytes",
    "Activated_Mfs", "Intermediate_Mfs", "Resident_Mfs", "Monocytes", "Cdc1", "Cdc2_a", 
    "Cdc2_b", "pDCs", "NKT_cells", "Cd8_T_cells", "Cycling_T_cells", "Memory_T_cells",
    "Naïve_CD4_T_cells", "Naïve_CD8_T_cells", "Tcr_gammadelta", "Th17_cells", "Tregs", "FDCs_MRCs",
    "FRC1", "FRC2", "Mesothelial", "SMCs", "Pericytes", "Cajal", 
    "Capillary_ECs", "LECS", "Vein_ECs", "Trophocytes", "Pdgfra_low", "Telocytes"
  )
)
)
dare_sum <- as.data.frame(dare_sum)
rownames(dare_sum) <- dare_sum$target

#Barplot for TNF receiver
bp_receiver <- ggplot(dare_sum, aes(target, prob, fill=target)) +  
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = as.vector(dare_sum[levels(dare_sum$target), 'colors']) ) +
  coord_flip()+
  theme_bw()+
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(angle = 65, hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 14)) +
  labs(x="Receiver", y="Probability", title = "TNF receiver")

#2)---source barplot---
dare_table <- subsetCommunication(cellchat.LS, slot.name = "netP", signaling = "TNF")
dare_sum2 <- dare_table %>% 
  group_by(source, pathway_name) %>% 
  summarise(prob = sum(prob))
df_label_colors2 <- data.frame(source=levels(gg[[1]]$data$labels),
                              colors=gg[[1]]$plot_env$color.use)
dare_sum2 <- left_join(df_label_colors2, dare_sum2)
dare_sum2$prob[is.na(dare_sum2$prob)] <- 0
dare_sum2$source <- factor(dare_sum2$source, levels = rev(
                                                        c("B_cells", "Cycling_B_cells", "Plasma_cells", "ILC2s", "ILC3s", "Granulocytes",
                                                        "Activated_Mfs", "Intermediate_Mfs", "Resident_Mfs", "Monocytes", "Cdc1", "Cdc2_a", 
                                                        "Cdc2_b", "pDCs", "NKT_cells", "Cd8_T_cells", "Cycling_T_cells", "Memory_T_cells",
                                                        "Naïve_CD4_T_cells", "Naïve_CD8_T_cells", "Tcr_gammadelta", "Th17_cells", "Tregs", "FDCs_MRCs",
                                                        "FRC1", "FRC2", "Mesothelial", "SMCs", "Pericytes", "Cajal", 
                                                        "Capillary_ECs", "LECS", "Vein_ECs", "Trophocytes", "Pdgfra_low", "Telocytes"
                                                         )
                                                        )
                          )
dare_sum2 <- as.data.frame(dare_sum2)
rownames(dare_sum2) <- dare_sum2$source

#Barplot for TNF source
bp_source <- ggplot(dare_sum2, aes(source, prob, fill=source)) +  
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = as.vector(dare_sum2[levels(dare_sum2$source), 'colors']) ) +
  coord_flip()+
  theme_bw()+
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_text(angle = 0), 
        plot.title = element_text(angle = 65, hjust = 0.5, vjust = 0.5), 
        text = element_text(size = 14)) +
  labs(y="Probability", title = "TNF source")

#3) heatmap of TNF activity
#calculation of activity scores for DARE sample

obj <- readRDS("../Step9_CellToCellCommunicationAnalysis/Stroma_Myeloid_Lymphoid_obj.RDS")
#----PROGENy model
net <- read.delim("utilityFunctions/progeny_top100_mouse.txt")

table(obj$sample) #keep dare-12w
Idents(obj) <- obj$sample
obj <- subset(obj, idents = "dare-12w")
Idents(obj) <- obj$annotation
table(obj$annotation)

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

#Get mean activity values
df <- t(as.matrix(obj@assays$pathwayswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(obj)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

#Dataframe used to plot mean scaled pathway activity values as heatmap
tnf_df <- df |> filter(source == "TNFa" & cluster != "Glial")
tnf_df$z_score <- (tnf_df$mean - mean(tnf_df$mean)) / sd(tnf_df$mean)
tnf_df$cluster <- factor(tnf_df$cluster, levels = rev(
  c("B_cells", "Cycling_B_cells", "Plasma_cells", "ILC2s", "ILC3s", "Granulocytes",
    "Activated_Mfs", "Intermediate_Mfs", "Resident_Mfs", "Monocytes", "Cdc1", "Cdc2_a", 
    "Cdc2_b", "pDCs", "NKT_cells", "Cd8_T_cells", "Cycling_T_cells", "Memory_T_cells",
    "Naïve_CD4_T_cells", "Naïve_CD8_T_cells", "Tcr_gammadelta", "Th17_cells", "Tregs", "FDCs_MRCs",
    "FRC1", "FRC2", "Mesothelial", "SMCs", "Pericytes", "Cajal", 
    "Capillary_ECs", "LECS", "Vein_ECs", "Trophocytes", "Pdgfra_low", "Telocytes"
  )
))

# Heatmap using ggplot
hmap <- ggplot(tnf_df, aes(x = source, y = cluster, fill = z_score)) +
  geom_tile(color = "black") +  
  scale_fill_gradient2(
    low = "darkblue",      
    mid = "white",   
    high = "red",    
    midpoint = 0,     
    limits = c(-1.5, 1.5),  
    oob = squish      
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),  
    panel.grid = element_blank(),  
    text = element_text(size = 14), 
    legend.position = "none",
    axis.text.x = element_blank(),
    plot.title = element_text(angle = 65, hjust = 0.5, vjust = 0.5)
  ) +
  ggtitle("TNF activity")

#save combined plot
combined_plot <- (hmap + bp_source + bp_receiver) + 
  plot_layout(widths = c(0.5, 1, 1))
combined_plot
ggsave("Fig4c.pdf", width = 14, height = 24, units = "cm", dpi = 600)
