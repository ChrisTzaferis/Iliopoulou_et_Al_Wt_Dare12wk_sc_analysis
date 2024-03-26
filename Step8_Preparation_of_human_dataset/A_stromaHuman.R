library(Seurat)
library(dplyr)

#Download the raw data from cellportal: https://singlecell.broadinstitute.org/single_cell/study/SCP1884/human-cd-atlas-study-between-colon-and-terminal-ileum#study-download 
#and store them in a folder named "stroma_with_raw/" iside the current directory 
#Download also the metadata file
str_data <- Read10X("stroma_with_raw/")
str_obj <- CreateSeuratObject(counts = str_data, project = "stroma", min.cells = 0, min.features = 0)
metadata <- read.delim("metadata.txt")

metadata_str <- metadata[which(metadata$organ__ontology_label == "ileum"), ]
metadata_str <- metadata_str[which(metadata_str$Celltype %in% c("Activated fibroblasts CCL19 ADAMADEC1", "Endothelial cells CA4 CD36", "Endothelial cells CD36", "Endothelial cells DARC",
                                                                "Endothelial cells LTC4S SEMA3G", "Fibroblasts ADAMDEC1", "Fibroblasts KCNN3 LY6H", "Fibroblasts NPY SLITRK6",
                                                                "Fibroblasts SFRP2 SLPI", "Fibroblasts SMOC2 PTGIS", "Glial cells", "Lymphatics",
                                                                "Myofibroblasts GREM1 GREM2", "Myofibroblasts HHIP NPNT", "Pericytes HIGD1B STEAP4", "Pericytes RERGL NTRK2"
)), ]
rownames(metadata_str) <- metadata_str$NAME

identical(rownames(str_obj@meta.data), rownames(metadata_str))#should be true
head(str_obj@meta.data)
head(metadata_str)

str_obj$NAME <- rownames(str_obj@meta.data)

#new metadata 
new_metadata <- left_join(str_obj@meta.data, metadata_str)
rownames(new_metadata) <- new_metadata$NAME
str_obj@meta.data <- new_metadata

saveRDS(str_obj, "TI_Stroma_Human.RDS")
