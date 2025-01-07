library(CellChat)
library(tidyverse)
library(ggplot2)
library(patchwork)

#Load objects
cellchat <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchatCombinedDARE_WT.RDS")
cellchat.NL <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("../Step9_CellToCellCommunicationAnalysis/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

#Keep signaling pathways where fibroblasts are considered as "sender" in DARE
dare_table <- subsetCommunication(cellchat.LS, slot.name = "net", sources.use = c("Pdgfra_low", "Trophocytes", "Telocytes"))
dare_table <- dare_table |> filter(pathway_name %in% c("CCL", "CXCL")) 
dare_sig_pathways <- unique(dare_table$interaction_name_2)
dare_table$sample <- "Dare-12w"
dare_table <- dare_table[, c(5, 8, 12)] #keep sample name, interaction name and probability of the interaction

#Keep signaling pathways where fibroblasts are considered as "sender" in WT
wt_table <- subsetCommunication(cellchat.NL, slot.name = "net", sources.use = c("Pdgfra_low", "Trophocytes", "Telocytes"))
wt_table <- wt_table |> filter(pathway_name %in% c("CCL", "CXCL")) 
wt_sig_pathways <- unique(wt_table$interaction_name_2)
wt_table$sample <- "Wt-12w"
wt_table <- wt_table[, c(5, 8, 12)]

#Aggregated probality when the three fibroblast clusters are considered as one group
wt_sum <- wt_table %>% 
  group_by(interaction_name_2, sample) %>% 
  summarise(prob = sum(prob))

dare_sum <- dare_table %>% 
  group_by(interaction_name_2, sample) %>% 
  summarise(prob = sum(prob))

#concatenate the two tables
final_table <- rbind(wt_sum, dare_sum)
head(final_table, n = 20)

# Common pathways and calculation of log2_prob_ratio
df_common <- final_table |>
  filter(interaction_name_2 %in% intersect(dare_sig_pathways, wt_sig_pathways) ) |>
  pivot_wider(names_from = sample, values_from = prob) %>%
  mutate(
    log2_prob_ratio = log2(`Dare-12w` / `Wt-12w`)
  ) %>%
  select(interaction_name_2, `Dare-12w`, `Wt-12w`, log2_prob_ratio)

#Dotplot
df_common <- df_common[order(df_common$log2_prob_ratio), ]
x_axis_colors <- c(rep("black", 7), "red4", "black", "red4", rep("black", 6), "red4", "red4", "black", "red4")

ggplot(df_common) + theme_bw() +
  geom_point(aes(x= "", y=reorder(interaction_name_2, log2_prob_ratio), size=5, color=log2_prob_ratio))+
  scale_colour_gradient2(low="blue", mid = "grey90", high = "red", midpoint = 0) +
  scale_size(range = c(9, 9), guide = "none") + 
  coord_flip() +
  theme(axis.text.x = element_text(face = "bold",  size = 16, color = x_axis_colors, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold", color = "black", size = 16, angle = 0),
        legend.title = element_text(face = "bold", color = "black", size = 16, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(size = 16, face = "bold.italic"),
        axis.ticks.y = element_blank()
  ) +
  ggtitle("Fibroblasts as senders") + 
  labs(x="", y="", color="log2(Prob. ratio)", size=" ")
ggsave("SFig11a.pdf", width = 25, height = 8, units = "cm", dpi = 600)
