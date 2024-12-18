library(Seurat)

#Lymphoid, Myeloid and Stroma merging
#------------load 3 objects------------------
myeloid <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")
DimPlot(myeloid)

lymphoid <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")
DimPlot(lymphoid)

stroma <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Stroma_annotated.RDS")
DimPlot(stroma)

colnames(myeloid@meta.data)
colnames(lymphoid@meta.data)
colnames(stroma@meta.data)

#merge objects
comb_obj <- merge(stroma, y = c(myeloid, lymphoid), add.cell.ids = c("str", "myel", "lym"), project = "dare12")
table(comb_obj$orig.ident)

#visualization of objects
comb_obj@assays$integrated@var.features <- union(VariableFeatures(stroma), c(VariableFeatures(lymphoid), VariableFeatures(myeloid)))
DefaultAssay(comb_obj) <- "integrated"
comb_obj <- ScaleData(comb_obj)
comb_obj <- RunPCA(comb_obj, npcs = 50)
ElbowPlot(comb_obj, ndims = 50)
comb_obj <- RunUMAP(comb_obj, dims = 1:30)
DimPlot(comb_obj, group.by = "annotation", label = T)

DefaultAssay(comb_obj) <- "RNA"
comb_obj <- ScaleData(comb_obj, features = VariableFeatures(comb_obj, assay = "integrated")) #keeping only mvgs at this point
saveRDS(comb_obj, "Stroma_Myeloid_Lymphoid_obj.RDS")

#Load obj and split to control and disease for cellchat analysis
init_obj <- readRDS("Stroma_Myeloid_Lymphoid_obj.RDS")
Idents(init_obj) <- init_obj$annotation

#Delete glial cells (only 1 cell in wt)
init_obj <- subset(init_obj, idents="Glial", invert=T)

#Split in two objects and save
Idents(init_obj) <- init_obj$sample
obj1 <- subset(init_obj, idents = "wt-12w") #control
obj2 <- subset(init_obj, idents = "dare-12w") #experimental
Idents(obj1) <- obj1$annotation
Idents(obj2) <- obj2$annotation
saveRDS(obj1, "WT12_cellchatInput.RDS")
saveRDS(obj2, "DARE12_cellchatInput.RDS")
