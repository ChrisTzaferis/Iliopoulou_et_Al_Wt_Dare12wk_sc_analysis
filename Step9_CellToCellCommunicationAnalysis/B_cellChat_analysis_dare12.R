library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

#DARE objcect containing Lymphoid, Myeloid, Stroma
curObj <- readRDS("DARE12_cellchatInput.RDS")

#Create cellchat object
curObj$annotation <- factor(curObj$annotation)
cellchat <- createCellChat(object = curObj@assays$RNA@data, meta = curObj@meta.data, group.by = "annotation")

#---DARE output---
# Create a CellChat object from a data matrix
# Set cell identities for the new CellChat object
# The cell groups used for CellChat analysis are  Activated_Mfs B_cells Cajal Capillary_ECs Cd8_T_cells Cdc1 
# Cdc2_a Cdc2_b Cycling_B_cells Cycling_T_cells FDCs_MRCs FRC1 FRC2 Granulocytes ILC2s ILC3s Intermediate_Mfs 
# LECS Memory_T_cells Mesothelial Monocytes Naïve_CD4_T_cells Naïve_CD8_T_cells NKT_cells pDCs Pdgfra_low 
# Pericytes Plasma_cells Resident_Mfs SMCs Tcr_gammadelta Telocytes Th17_cells Tregs Trophocytes Vein_ECs

cellchat <- addMeta(cellchat, meta = curObj@meta.data)
cellchat <- setIdent(cellchat, ident.use = "annotation") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

#---DARE output---
# [1] "Activated_Mfs"     "B_cells"           "Cajal"             "Capillary_ECs"     "Cd8_T_cells"       "Cdc1"              "Cdc2_a"           
# [8] "Cdc2_b"            "Cycling_B_cells"   "Cycling_T_cells"   "FDCs_MRCs"         "FRC1"              "FRC2"              "Granulocytes"     
# [15] "ILC2s"             "ILC3s"             "Intermediate_Mfs"  "LECS"              "Memory_T_cells"    "Mesothelial"       "Monocytes"        
# [22] "Naïve_CD4_T_cells" "Naïve_CD8_T_cells" "NKT_cells"         "pDCs"              "Pdgfra_low"        "Pericytes"         "Plasma_cells"     
# [29] "Resident_Mfs"      "SMCs"              "Tcr_gammadelta"    "Telocytes"         "Th17_cells"        "Tregs"             "Trophocytes"      
# [36] "Vein_ECs"   

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

CellChatDB <- CellChatDB.mouse #if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

future::plan("multisession", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

library(NMF)
k1 <- selectK(cellchat, pattern = "outgoing") #
ggsave(filename="k_outgoing_DARE12.png", plot = k1, width=12, height=6, units = "cm")

k2 <- selectK(cellchat, pattern = "incoming") #
ggsave(filename="k_incoming_DARE12.png", plot = k2, width=12, height=6, units = "cm")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional", do.parallel=F)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural", do.parallel=F)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
saveRDS(cellchat, "cellchat_DARE12_final.RDS")
