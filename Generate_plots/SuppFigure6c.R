library(dplyr)
library(viridis)
library(ggplot2)

#file produced from metascape
meta_table <- read.csv("../Metascape_files/sfig6c_GO_AllLists.csv")

#selected terms
selection <- c("inflammatory response", "positive regulation of cell migration", 
               "cell chemotaxis", "positive regulation of response to external stimulus", 
               "response to interferon-beta", "response to tumor necrosis factor", 
               "regulation of endothelial cell proliferation", "negative regulation of lipid localization", 
               "response to copper ion", "Non-integrin membrane-ECM interactions", "regulation of collagen biosynthetic process", 
               "Collagen degradation", "inflammatory response to wounding", "response to oxygen levels", 
               "antigen processing and presentation of exogenous peptide antigen", "supramolecular fiber organization", 
               "leukocyte activation", "cell-cell adhesion", "regulation of neuron projection development", 
               "positive regulation of cell development")

#keep selected terms and enrichment values
meta_df <- meta_table |> filter(Description %in% selection) |> 
  select(Description, Enrichment, GeneList, Log.q.value.) |> 
  mutate(LogP = -Log.q.value.)

#Order of plotting
meta_df$Description <- factor(meta_df$Description, levels = selection)
meta_df$GeneList <- gsub(pattern = "_UP", replacement = "", x = meta_df$GeneList)
meta_df$GeneList <- factor(meta_df$GeneList, levels = c("Pdgfra_low", "Telocytes", "Trophocytes", "S1", "S2", "S3"))                                          

#Dotplot
ggplot(meta_df, aes(x = GeneList, y= Description, size=Enrichment, fill= LogP)) +
  geom_point(shape = 21, colour="black", stroke=0.5) +
  scale_size_continuous(range = c(1.5, 8))+
  scale_fill_viridis(option="F") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), axis.title = element_blank()) +
  labs(title = "", x = "Clusters", y = "Term", fill="-log10(P adj.)")+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="grey")))
ggsave(filename = "SFig6c.pdf", width = 18, height = 15, units = "cm", dpi = 600)