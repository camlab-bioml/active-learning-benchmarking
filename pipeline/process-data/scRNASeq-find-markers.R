suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DESeq2)
  library(tidyverse)
  library(annotables)
  library(scater)
  library(scran)
  library(ComplexHeatmap)
})
source("pipeline/whatsthatcell-helpers.R")

sce <- readRDS(snakemake@input$sce)

fm <- findMarkers(sce, sce$CellType)

top_markers_per_cluster <- lapply(fm, function(f) {
  rownames(head(f, n = 20))
})

top_markers <- unique(unlist(top_markers_per_cluster))

expr_mat <- sce[top_markers,] %>% 
  logcounts() %>% 
  t() %>% 
  scale() 

thresh <- 2
expr_mat[expr_mat > thresh] <- thresh
expr_mat[expr_mat < -thresh] <- -thresh

dfc <- expr_mat %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(cluster = sce$CellType) %>% 
  gather(gene, expression, -cluster)

dfc2 <- group_by(dfc, cluster, gene) %>% 
  summarise(mean_expression = mean(expression), 
            pct_expressing = 100 * mean(expression > 0)) %>% 
  ungroup()


hm_mat <- dfc2 |> 
  select(-pct_expressing) |> 
  pivot_wider(names_from = "gene", values_from = "mean_expression") |> 
  as.data.frame() |> 
  column_to_rownames("cluster") |> 
  as.matrix()
  
hm <- Heatmap(hm_mat, name = "lfc")
hm <- draw(hm)

pdf(snakemake@output$heatmap, height = 10, width = 16)
  hm
dev.off()

cell_line_order <- rownames(hm_mat)[row_order(hm)]
gene_order <- colnames(hm_mat)[column_order(hm)]

pdf(snakemake@output$lfc_point, height = 10, width = 16)
  dfc2 |> 
    mutate(cluster = factor(cluster, levels = cell_line_order),
           gene = factor(gene, levels = gene_order)) |> 
  ggplot(aes(y = cluster, x = gene)) +
    geom_point(aes(colour = mean_expression, size = pct_expressing)) +
    scale_colour_viridis_c(name = "Expression") +
    scale_size(name = "%  cells expressing") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top") +
    labs(x = "Gene", y = "Cluster")
dev.off()



