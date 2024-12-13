library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library("ggsci")
library(scales)
library(grid)
library(viridis)

#Enriched terms in metascape
all_terms <- read.delim("../Metascape_files/fig3d_GO_AllLists.csv")

#Selected terms to be plotted
sel_terms <- c("antigen processing and presentation", "response to type II interferon", "response to interleukin-1", 
                              "TNF signaling pathway - Mus musculus (house mouse)", "neutrophil migration", "NOD-like receptor signaling pathway - Mus musculus (house mouse)", 
                              "positive regulation of cell migration", "regulation of leukocyte migration", "inflammatory response", "response to lipopolysaccharide", 
                              "cellular response to lipid", "regulation of chemotaxis", "positive regulation of ERK1 and ERK2 cascade", "response to wounding", 
                              "Extracellular matrix organization", "blood vessel morphogenesis", "vasculature development", "tube morphogenesis", "blood circulation")

all_terms_filt <- all_terms[which(all_terms$Description %in% sel_terms), ] #keep selected terms

all_terms_filt$Description <- factor(all_terms_filt$Description, levels = sel_terms)
all_terms_filt$Count <- all_terms_filt$X.GeneInGOAndHitList #gene count (gene overlap between input list and term)
all_terms_filt$Cluster <- gsub(all_terms_filt$GeneList, pattern = "_UP", replacement = "", perl = T) # fix cluster names
all_terms_filt$Cluster <- factor(all_terms_filt$Cluster, levels = c("SMCs", "Pdgfra_low", "Trophocytes", "Telocytes","Capillary_ECs",
                                                                    "LECS", "Vein_ECs","Pericytes","FDCs_MRCs","Mesothelial",
                                                                    "FRC1", "FRC2", "Cajal"))

#Enrichment values greater than 6 will be set to 6 and will have the same color in the dotplot
for(i in 1:length(all_terms_filt$Enrichment))
{
  if(all_terms_filt$Enrichment[i] > 6)
  {
    all_terms_filt$Enrichment[i] <- 6
  }
}

#Print and save dotplot
ggplot(all_terms_filt) + theme_bw() +
  geom_point(aes(x=Cluster, y=Description, size=Count , color=Enrichment))+
  scale_size(range = c(2, 12), breaks = c(5, 15, 30)) +    
  scale_color_viridis( option = "C") +
  facet_grid(~Cluster, scales = "free_x") +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10, angle = 30, vjust = 0.8, hjust = 0.9),
        axis.text.y = element_text(face = "bold", color = "black", size = 14, angle = 0),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="white", fill="white", size=1, linetype="solid")) +
  labs(color="Enrichment", size="Gene \ncount")
ggsave("Fig3e.pdf", device = "pdf", width = 14, height = 7, units = "in", dpi = 600)
