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

#Keep signaling pathways where fibroblasts are considered as "receiver" in DARE
dare_table <- subsetCommunication(cellchat.LS, slot.name = "netP", targets.use = c("Pdgfra_low", "Trophocytes", "Telocytes"))
dare_sig_pathways <- unique(dare_table$pathway_name)
dare_table$sample <- "Dare-12w"
dare_table <- dare_table[, c(3, 4, 6)]

#Keep signaling pathways where fibroblasts are considered as "receiver" in WT
wt_table <- subsetCommunication(cellchat.NL, slot.name = "netP", targets.use = c("Pdgfra_low", "Trophocytes", "Telocytes"))
wt_sig_pathways <- unique(wt_table$pathway_name)
wt_table$sample <- "Wt-12w"
wt_table <- wt_table[, c(3, 4, 6)]

#Aggregated probality when the three fibroblast clusters are considered as one group
wt_sum <- wt_table %>% 
  group_by(pathway_name, sample) %>% 
  summarise(prob = sum(prob))

dare_sum <- dare_table %>% 
  group_by(pathway_name, sample) %>% 
  summarise(prob = sum(prob))

#concatenate the two tables
final_table <- rbind(wt_sum, dare_sum)
head(final_table, n = 20)

#unique pathways to Dare
df_dare12w <- final_table |>
  filter(pathway_name %in% setdiff(dare_sig_pathways, wt_sig_pathways))

#unique pathways to Wt
df_wt12w <- final_table |>
  filter(pathway_name %in% setdiff(wt_sig_pathways, dare_sig_pathways))

# Common pathways and calculation of log2_prob_ratio
df_common <- final_table |>
  filter(pathway_name %in% intersect(dare_sig_pathways, wt_sig_pathways) ) |>
  pivot_wider(names_from = sample, values_from = prob) %>%
  mutate(
    log2_prob_ratio = log2(`Dare-12w` / `Wt-12w`)
  ) %>%
  select(pathway_name, `Dare-12w`, `Wt-12w`, log2_prob_ratio)

#Dotplots
#wt dotplot
df_wt12w$prob[df_wt12w$prob > 0.033] <- 0.034
p1 <- ggplot(df_wt12w) + theme_bw() +
  geom_point(aes(x= sample, y=reorder(pathway_name, prob), size=5, color=prob))+
  scale_colour_gradient2(low="white", mid = "grey90", high = "darkblue") +
  scale_size(range = c(9, 9), guide = "none") + 
  coord_flip() + 
  theme(axis.text.x = element_text(face = "bold", color = "darkblue", size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(face = "bold", color = "black", size = 18, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(size = 16)
  ) +
  ggtitle(expression(bolditalic('Unique in Tnf'^'+/+')))+
  labs(x="", y="", color="Probability", size=" ")
p1

#dare dotplot
p2 <- ggplot(df_dare12w) + theme_bw() +
  geom_point(aes(x= sample, y=reorder(pathway_name, prob), size=5, color=prob))+
  scale_colour_gradient2(low="white", mid = "grey90", high = "red") +
  scale_size(range = c(9, 9), guide = "none") + 
  coord_flip() + 
  theme(axis.text.x = element_text(face = "bold", color = "red", size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(face = "bold", color = "black", size = 18, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(size = 16)
  ) +
  ggtitle(expression(bolditalic('Unique in Tnf'^'DARE')))+
  labs(x="", y="", color="Probability", size=" ")
p2

#Dotplot of shared signaling pathways
df_common <- df_common[order(df_common$log2_prob_ratio), ]

x_axis_colors <- ifelse(
  df_common$log2_prob_ratio < -0.15680827, "blue",
  ifelse(df_common$log2_prob_ratio < 0.17807123, "grey", "red")
)

p3 <- ggplot(df_common) + theme_bw() +
  geom_point(aes(x= "", y=reorder(pathway_name, log2_prob_ratio), size=5, color=log2_prob_ratio))+
  scale_colour_gradient2(low="blue", mid = "grey90", high = "red", midpoint = 0) +
  scale_size(range = c(9, 9), guide = "none") + 
  coord_flip() +
  theme(axis.text.x = element_text(face = "bold", color = x_axis_colors, size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold", color = "black", size = 16, angle = 0),
        legend.title = element_text(face = "bold", color = "black", size = 16, angle = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 12),
        axis.title.x = element_text(face = "bold", color = "black", size = 12),
        plot.title = element_text(size = 16, face = "bold.italic"),
        axis.ticks.y = element_blank()
  ) +
  ggtitle("Fibroblasts as receivers") + 
  labs(x="", y="", color="log2(Prob. ratio)", size=" ")
p3

#Combine all dotplots
combined_plot <- p3 / (p1 + p2 + plot_spacer()) +
  plot_layout(heights = c(1, 0.7), widths = c(1, 0.2))
ggsave(plot = combined_plot, "SFig7d.pdf", width = 52, height = 15, units = "cm", dpi = 600)
