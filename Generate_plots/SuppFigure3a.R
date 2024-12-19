library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

#---load lymphoid object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
Idents(obj) <- obj$annotation

#---Selected markers---
sel_genes <- c("S100a8", "Lrg1", "Chil1", "Cxcr2",
               "Chil3", "Fn1", "Ly6c2", "Ccr2",
               "Thbs1", "Fos", "Apoe", "C1qa",
               "Cx3cr1", "Acp5", "Saa3", "Lyz1",
               "Gpnmb", "Ccl5", "Xcr1", "Tlr3",
               "Ctla2b", "Clec9a", "Plet1", "H2-DMb2",
               "Rtn1", "Irf4", "Syn3", "Tcf7", "Dpp10",
               "Dscam", "Bst2", "Grm8", "Ccr9", "Siglech"
               )

#Plot and save
DotPlot(obj, features = sel_genes, dot.scale = 13, group.by = "annotation", assay = "RNA", scale = T) +
  theme_classic() +
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="plasma") +
  theme(text = element_text(size = 25), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )+
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="grey")))
ggsave(filename = "SFig3a.pdf", device = "pdf", width = 49, height = 13, units = "cm", dpi = 600)
