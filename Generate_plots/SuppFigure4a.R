library(Seurat)
library(dplyr)
library(ggplot2)

#---load lymphoid object---
obj <- readRDS("../Step7_Annotation_of_clusters/WT12wk_DARE12Wk_Myeloid_annotated.RDS")

#---DEGs calculation---
obj$cluster_sample <- paste0(obj$annotation, "_", obj$sample)

#clusters
cl_names <- as.vector(unique(obj$annotation))
#samples
sam_names <- unique(obj$sample)
#check cluster_sample column and exclude clusters: 1) containing less than 3 cells, 2) present only in one 
temp_df <- as.data.frame(table(obj$cluster_sample))
temp_df$Var1 <- gsub(pattern = "_dare-12w|_wt-12w", replacement = "", x = temp_df$Var1)
temp_df <- temp_df |> group_by(Var1) |> mutate(n=n()) #n=2 cluster is found in both conditions, n=1 cluster is present only in one condition
excl_clusters <- as.vector(temp_df[which(temp_df$Freq < 3 | temp_df$n != 2), "Var1"])

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

#--------------------Up and down regulated genes--------------------------------
up <- df[which(df$avg_logFC > 0), ]
down <- df[which(df$avg_logFC < 0), ]

#------------------------As a table---------------------------------------------
total_up <- up %>% count(cluster, name = "UP")
total_down <- down %>% count(cluster, name = "DOWN")

total_degs <- left_join(total_up, total_down)
total_degs$TOTAL <- total_degs$UP+ total_degs$DOWN
total_degs$TYPE <- "Myeloid"

#-----------------------for R barplots------------------------------------------
#modify the dataframe and insert to ggplot
total_up_down_bar <- rbind(
  total_degs |> select(cluster, UP) |> mutate(state="up", n=UP) |> select(-UP), 
  total_degs |> select(cluster, DOWN) |> mutate(state="down", n=DOWN) |> select(-DOWN)
)

#reorder clusters based on total number of DEGs
total_up_down_bar$cluster <- factor(total_up_down_bar$cluster, levels = total_degs$cluster[order(total_degs$TOTAL, decreasing = F)])
total_up_down_bar$state <- factor(total_up_down_bar$state, levels = c("down", "up"))

ggplot(total_up_down_bar, aes(fill=state, y=n, x=cluster)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c("royalblue", "red3"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), text = element_text(size = 29), plot.title = element_text(hjust = 0.5)) +
  labs(y="n of DEGs", x="", title = expression(italic(Tnf^{DARE}) ~ "vs" ~ italic(Tnf^{"+/+"})) )
ggsave("SFig4a.pdf", width = 18, height = 14, units = "in", dpi = 600)

