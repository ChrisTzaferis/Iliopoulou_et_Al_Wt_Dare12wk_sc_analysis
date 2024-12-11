library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

# meta_table: metascape file GO_AllLists.csv
# selectedTerms: selected terms in a vector. The names should be included in the "Description" column in the GO_AllLists.csv file
# selectedClusters: selected clusters in a vector. The names should be included in the "GeneList" column in the GO_AllLists.csv file
# maxCol= maximum value for color itensity, higher values will be set to maxCol.
# pdfName, pdfWidth, pdfHeight: export options
# clusterRows, clusterCols: enable clustering of rows/columns

metascapeToHeatmap <- function(meta_table, selectedTerms, selectedClusters, maxCol=100, pdfName, pdfWidth=12, pdfHeight=12, clusterRows=F, clusterCols=F)
{
  #keep selected terms and enrichment values
  meta_df <- meta_table |> 
    filter(Description %in% selectedTerms & GeneList %in% selectedClusters) |> 
    select(Description, Enrichment, GeneList) |>
    mutate(Enrichment = ifelse(Enrichment > maxCol, maxCol, Enrichment))
  
  #convert to matrix
  heatmap_table <- meta_df %>%
    pivot_wider(id_cols = 'Description', names_from = 'GeneList',
                values_from = 'Enrichment') %>%
    tibble::column_to_rownames('Description') %>%
    as.matrix()
  
  #replace NA with 0
  heatmap_table[is.na(heatmap_table)] <- 0
  
  #print heatmap
  pdf(file = pdfName, height = pdfHeight, width = pdfWidth)
  pheatmap(heatmap_table[selectedTerms, selectedClusters], 
           legend = T,
           cluster_rows = clusterRows, 
           cluster_cols = clusterCols, 
           color = colorRampPalette(c("white", "lightpink","red", "red3", "red4"))(30))
  dev.off()
}