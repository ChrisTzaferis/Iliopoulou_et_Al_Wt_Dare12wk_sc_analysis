library(tidyr)
library(dplyr)
library(Seurat)
library(scGSVA)  
library(pheatmap)
library(ggplot2)
library(viridis)

# seu_obj: a seurat object
# geneSetType: KEGG/REACTOME/HALLMARK/GO/BP/CC/MF
# calcMethod: ssgsea/UCell
performSCGSVA <- function(seu_obj, geneSetType, calcMethod)
{
  set.seed(123)
  
  #get genesets
  gs <- buildMSIGDB(species="mouse",keytype = "SYMBOL",anntype = geneSetType)
  
  #calculate enrichment scores
  res <- scgsva(seu_obj, annot = gs, method = calcMethod, assay = "RNA", slot = "counts", verbose = T)
  
  #store scores in object
  seu_obj[["ssGSEA"]] <- CreateAssayObject(counts = t(as.data.frame(res)))
  
  #return object containing the GSEA assay
  return(seu_obj)
}

# seu_obj: a seurat object 
# groupCluster: a metadata column used to group cells before calculations
# assayName: an assay contatining the results from scGSVA
# identsToExclude: remove clusters from visualization
# showAll: plot all terms
# nToPlot: plot the top 'nToPlot' terms per cluster ranked by mean activity
# filenameOut: save heatmap under this filename
# imgWidth, imgHeight: dimensions of the plot

plotTopSD <- function(seu_obj, groupCluster="cluster_sample", assayName="ssGSEA", identsToExclude=NULL, showAll=F, nToPlot=5, filenameOut, imgWidth=30, imgHeight=40)
{
  
  Idents(seu_obj) <- groupCluster
  DefaultAssay(seu_obj) <- assayName
  if(!is.null(identsToExclude) )
  {
    if( identsToExclude %in% unique(Idents(seu_obj)) )
      seu_obj <- subset(x = seu_obj, idents = identsToExclude, invert=T)
  }
  
  seu_obj <- ScaleData(seu_obj, assay = assayName)
  
  # Extract GSEA scores from object as a long dataframe
  df <- t(as.matrix(seu_obj@assays[[assayName]]@scale.data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(seu_obj)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score))
  
  if(showAll == T)
  {
    terms <- unique(df$source)
  }else
  {
    terms <- df |> 
      group_by(cluster) |> 
      top_n(n = nToPlot, wt = mean) |> 
      pull(source) |>
      unique() |>
      as.character() |>
      as.vector()
  }
  
  # Subset long data frame to top terms and transform to wide matrix
  top_acts_mat <- df %>%
    filter(source %in% terms) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    tibble::column_to_rownames('cluster') %>%
    as.matrix()
  
  # create z-scores across clusters for all different terms
  scaled_mat <- scale(top_acts_mat)
  
  #min, max values in heatmap
  minVal <- min(scaled_mat)
  maxVal <- max(scaled_mat)
  
  #create color scheme
  palette_length = 100
  my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(minVal, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, maxVal, length.out=floor(palette_length/2)))
  
  #utilize pheatmap
  ph <- pheatmap(mat = t(scaled_mat), 
                 color = my_color, 
                 breaks = my_breaks, 
                 cluster_rows = T, 
                 cluster_cols = T,
                 display_numbers = F,
                 show_rownames = T,
                 border_color = "black"
  )
  ph
  
  #plot heatmap using ggplot2
  top_acts_mat2 <- scaled_mat %>%
    as.data.frame() %>%
    mutate(cluster = rownames(scaled_mat)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score")
  
  ordered_terms <- colnames(scaled_mat)[ph$tree_row$order]
  ordered_clusters <- rownames(scaled_mat)[ph$tree_col$order]
  top_acts_mat2$cluster <- factor(top_acts_mat2$cluster, levels = ordered_clusters)
  top_acts_mat2$source <- factor(top_acts_mat2$source, levels = rev(ordered_terms))
  
  ggplot(top_acts_mat2, 
         aes(cluster, source, fill= score)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          ) + labs(fill="z-score")
  
  ggsave(filename = filenameOut, width = imgWidth, height = imgHeight, units = "cm")
} 

