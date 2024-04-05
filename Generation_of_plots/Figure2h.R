library(Seurat)
library(dittoSeq)
library(Scillus)
library(UCell)
library(ggplot2)

#Load object
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
final_order <- c("Granulocytes","Monocytes", "Intermediate_Mfs", "Resident_Mfs",
                 "Activated_Mfs","Cdc1","Cdc2_a", "Cdc2_b","pDCs" )
obj$annotation <- factor(obj$annotation, levels = final_order)
obj$sample <- factor(obj$sample, levels = c("dare-12w", "wt-12w"))

#keep macrophages 
obj <- subset(obj, idents = c("Monocytes", "Intermediate_Mfs", "Activated_Mfs", "Resident_Mfs"))

#Signatures from GIMATs publication
#Human genes were converted to mouse genes using gProfiler
res_sig <- c("Folr2", "Fuca1", "Gpnmb", "Apoc1", "Lgmn", "Stab1", "Csf1r", "Ms4a4a", "Slc40a1", 
             "Slco2b1", "Mafb", "Mrc1", "Dab2", "Dnase1l3", "Rgs2", "Jaml", "Vsig4", "Cpvl", "Clec10a", 
             "Rnase6", "C1qa", "C1qb", "C1qc", "Cd14", "Cd68", "Pld3")

act_sig <- c("Cxcl2", "Cxcl3", "Cd83", "Sod2", "Txn1", "Ccl3", "Ccl4", "Nfkbia", "Cd44", "Pkm", "S100a6", "Il1b", 
             "Plaur", "Ier3", "G0s2", "S100a9", "S100a8", "Ccl20", "Ido1", "Inhba", "Il1rn", "Tnf", "Tnfaip6", 
             "Cxcl10", "Cxcl9", "Ninj1", "Il6", "Capg", "Aqp9", "Mmp9", "Ldha", "Kynu", "Il23a", "Wtap", "Ptgs2", 
             "Traf1", "Plek", "Gk", "Il1a", "Nfkb1", "Pim3", "Cflar", "Il4i1", "Birc3", "Snx10", "Cd40", "Icam1", 
             "Tspo", "Pde4dip", "Tnfaip2", "Gbp5", "Dusp2", "Cdkn1a", "Irf1")

#calculate signature scores
obj <- AddModuleScore_UCell(obj, features= list(Resident_Mfs_signature = res_sig, ProInfl_Mfs_signature = act_sig), assay = "RNA")

obj$annotation <- factor(obj$annotation, levels = c("Monocytes", "Intermediate_Mfs", "Resident_Mfs", "Activated_Mfs"))
Idents(obj) <- obj$annotation
p1 <- VlnPlot(obj, features =c("Resident_Mfs_signature_UCell") , pt.size = 0, ncol = 1, split.by = "sample", cols=c("red3", "grey")) + 
  geom_hline(yintercept = 0.3, linetype="dashed", color = "black", linewidth = 1) +
  ggtitle("Human resident Mfs signature (HC/CD)") +
  theme(axis.title.x = element_blank())

p2 <- VlnPlot(obj, features =c("ProInfl_Mfs_signature_UCell") , pt.size = 0, ncol = 1, split.by = "sample", cols=c("red3", "grey")) + 
  geom_hline(yintercept = 0.35, linetype="dashed", color = "black", linewidth = 1) +
  ggtitle("Human pro-inflammatory Mfs signature (HC/CD)")  +
  theme(axis.title.x = element_blank())

cowplot::plot_grid(p1, p2,
                   ncol = 1, align = "v")
ggsave("Fig2h.pdf", device = "pdf", width = 18, height = 25, units = "cm", dpi = 600)
