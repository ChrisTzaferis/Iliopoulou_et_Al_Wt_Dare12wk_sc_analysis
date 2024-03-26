library(Seurat)
library(dplyr)
library(nichenetr)

convertObject <- function(obj_path, proj_name, out_name)
{
  #---Convert human genes to mouse genes---
  current_obj <- readRDS(obj_path)
  genes_human <- rownames(current_obj)
  
  #from nichenetR
  genes_mouse <- genes_human %>% convert_human_to_mouse_symbols()
  df_mouse_human <- data.frame(genes_human, genes_mouse)
  df_mouse_human_final <- na.omit(df_mouse_human)
  length(intersect(df_mouse_human_final$genes_human, genes_human))
  colnames(df_mouse_human_final) <- c("genes_human", "genes_mouse")
  #identical(df_mouse_human$genes_human, rownames(current_obj))
  
  #use the genes from conversion
  obj_new <- subset(current_obj, features = df_mouse_human_final$genes_human)
  identical(df_mouse_human_final$genes_human, rownames(obj_new@assays$RNA@counts))
  max(obj_new@assays$RNA@counts["CCL19", ])
  
  rownames(obj_new@assays$RNA@counts) <- df_mouse_human_final$genes_mouse
  max(obj_new@assays$RNA@counts["Ccl19", ])
  
  final_new_obj <- CreateSeuratObject(counts = obj_new@assays$RNA@counts, assay = "RNA", meta.data = obj_new@meta.data, project = proj_name)
  saveRDS(final_new_obj, paste0(out_name, ".RDS"))
}

convertObject("./TI_Stroma_Human.RDS", "TI_Stroma", "TI_Stroma_Human_with_mouse_genes")
convertObject("./TI_Lymphoid_Human.RDS", "TI_Lymphoid", "TI_Lymphoid_Human_with_mouse_genes")
convertObject("./TI_Myeloid_Human.RDS", "TI_Myeloid", "TI_Myeloid_Human_with_mouse_genes")
