library(Seurat)
library(dplyr)
library(ggplot2)

#---load lymphoid object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Lymphoid_annotated.RDS")

#---DEGs calculation---
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)

#clusters
cl_names <- as.vector(unique(obj$annotation))
#samples
sam_names <- unique(obj$sample)
#check cluster_sample column and exclude small size clusters
temp_df <- as.data.frame(table(obj$cluster_sample))
temp_df$Var1 <- gsub(pattern = "_dare-12w|_wt-12w", replacement = "", x = temp_df$Var1)
excl_clusters <- as.vector(temp_df[which(temp_df$Freq < 3), "Var1"])

#clusters to be used for intra-cluster DEA
cl_names <- setdiff(cl_names, excl_clusters)

#define experimental conditions, control conditions as paired vectors
exp_vec <- paste0(cl_names, "_", sam_names[2])
ctr_vec <- paste0(cl_names, "_", sam_names[1])

#set "cluster_sample" column as Idents of the object
Idents(obj) <- obj$cluster_sample
DimPlot(obj)

#start intra-cluster DEA 
df <- data.frame()
for(i in 1:length(exp_vec))
{
  print(paste0(exp_vec[i], "_VS_", ctr_vec[i]))
  cur_df <- FindMarkers(obj, assay = "RNA", ident.1 = exp_vec[i], ident.2 = ctr_vec[i], logfc.threshold = 0.25, min.pct = 0.1, base = exp(1))
  cur_df$gene <- rownames(cur_df)
  cur_df <- cur_df[which(cur_df$p_val < 0.01), ]
  cur_df$comparison <- paste0(exp_vec[i], "_VS_", ctr_vec[i])
  cur_df$cluster <- gsub(pattern = "_dare-12w", replacement = "", x = exp_vec[i])
  df <- rbind(df, cur_df)
}

#---Barplot plotting---
up <- df[which(df$avg_logFC > 0), ]
down <- df[which(df$avg_logFC < 0), ]

#------------------------As a table---------------------------------------------
total_up <- up %>% count(cluster, name = "UP")
total_down <- down %>% count(cluster, name = "DOWN")

total_degs <- left_join(total_up, total_down)
total_degs$TOTAL <- total_degs$UP+ total_degs$DOWN
total_degs$TYPE <- "Lymphoid"

#-----------------------for R barplots------------------------------------------
total_up_bar <- up %>% count(cluster)
total_up_bar$state <- "up"
total_up_bar$type <- "lymphoid"
total_down_bar <- down %>% count(cluster)
total_down_bar$state <- "down"
total_down_bar$type <- "lymphoid"

total_up_down_bar <- rbind(total_up_bar, total_down_bar)
total_up_down_bar$cluster <- factor(total_up_down_bar$cluster, levels = total_up_bar$cluster[order(total_degs$TOTAL, decreasing = F)]) #ordered by total
total_up_down_bar$state <- factor(total_up_down_bar$state, levels = c("down", "up"))

ggplot(total_up_down_bar, aes(fill=state, y=n, x=cluster)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c("royalblue", "red3"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0), text = element_text(size = 20)) +
  coord_flip()+
  labs(y="n of DEGs", x="", title = "DARE vs Wt")
ggsave("Fig1g.pdf", width = 16, height = 15, units = "in", dpi = 600)

