library(Seurat)
library(dplyr)

imm_data <- Read10X("immune_with_raw/")
imm_obj_raw <- CreateSeuratObject(counts = imm_data, project = "immune", min.cells = 0, min.features = 0)

#ileum cells metadata
metadata <- read.delim("metadata.txt")
metadata_imm <- metadata[which(metadata$organ__ontology_label == "ileum"), ]

#lymphoid metadata
metadata_lymph <- metadata_imm[which(metadata_imm$Celltype %in% c(
  "B cells", "B cells AICDA LRMP", "IELs ID3 ENTPD1",
  "ILCs", "NK cells KLRF1 CD3G-", "Plasma cells", 
  "T cells CD4 FOSB", "T cells CD4 IL17A", "T cells CD8",
  "T cells CD8 KLRG1", "T cells Naive CD4", "T cells OGT",
  "Tregs"
)), ]
table(metadata_lymph$Celltype)
rownames(metadata_lymph) <- metadata_lymph$NAME

#myeloid metadata
metadata_myel <- metadata_imm[
  which(metadata_imm$Celltype %in% c(
    "Cycling cells", "DC1", "DC2 CD1D",
    "DC2 CD1D-", "Macrophages", "Macrophages CCL3 CCL4", 
    "Macrophages CXCL9 CXCL10", "Macrophages LYVE1", "Macrophages PLA2G2D",
    "Mature DCs", "Monocytes CHI3L1 CYP27A1", "Monocytes HBB",
    "Monocytes S100A8 S100A9", "Neutrophils S100A8 S100A9"
  )), ]
table(metadata_myel$Celltype)
rownames(metadata_myel) <- metadata_myel$NAME

length(which(rownames(imm_obj_raw@meta.data) %in% rownames(metadata_lymph))) #143552
length(which(rownames(imm_obj_raw@meta.data) %in% rownames(metadata_myel))) #48980

#lymph object
obj_lymph <- subset( imm_obj_raw, cells = rownames(metadata_lymph) )
new_lymph_meta <- obj_lymph@meta.data
new_lymph_meta$NAME <- rownames(new_lymph_meta)
new_lymph_meta <- left_join(new_lymph_meta, metadata_lymph)
rownames(new_lymph_meta) <- new_lymph_meta$NAME
obj_lymph@meta.data <- new_lymph_meta
head(obj_lymph@meta.data)

#myeloid object
obj_myel <- subset( imm_obj_raw, cells = rownames(metadata_myel) )
new_myel_meta <- obj_myel@meta.data
new_myel_meta$NAME <- rownames(new_myel_meta)
new_myel_meta <- left_join(new_myel_meta, metadata_myel)
rownames(new_myel_meta) <- new_myel_meta$NAME
obj_myel@meta.data <- new_myel_meta
head(obj_myel@meta.data)

saveRDS(obj_lymph, file = "TI_Lymphoid_Human.RDS")
saveRDS(obj_myel, file = "TI_Myeloid_Human.RDS")
