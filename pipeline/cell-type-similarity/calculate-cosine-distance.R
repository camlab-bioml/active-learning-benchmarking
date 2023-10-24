suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
  library(data.table)
  library(RPortfolioSimilarity)
  library(ComplexHeatmap)
})
source("pipeline/whatsthatcell-helpers.R")


sce <- readRDS(snakemake@input$sce)

if(!is.null(snakemake@params$pca)){
  df_PCA <- reducedDim(sce, "PCA")
  var_explained <- tibble(
    PC = colnames(df_PCA),
    `Proportion of Variance` = attr(reducedDim(sce, "PCA"), "percentVar")
  )

  df_PCA <- df_PCA[,1:20] |>
    as.data.frame() |>
    mutate(cell_type = sce$CellType)
  var_explained <- var_explained[1:20, ]
}else{
  if(is.null(sce$CellType)){
    sce$CellType <- sce$cell_type
  }

  df_expression <- load_scs(sce) |> 
    select(-X1)

  # Remove genes with 0 expression
  expression <- df_expression[, colSums(df_expression) > 0]
  # Run PCA
  df_PCA <- as.matrix(expression) |> 
    prcomp(center = TRUE, scale. = TRUE)

  # Get variance explained
  var_explained <- summary(df_PCA)$importance |> 
    t() |>
    as.data.frame() |> 
    rownames_to_column("PC") |> 
    select(PC, "Proportion of Variance") 
  var_explained <- var_explained[1:20,]

  # Get PCA embedding
  df_PCA <- df_PCA$x[,1:20] |> 
    as.data.frame() |> 
    mutate(cell_type = sce$CellType)
}

celltypes <- unique(df_PCA$cell_type)

cell_type_sims <- lapply(celltypes, function(cell_type_i){
  sim <- lapply(celltypes, function(cell_type_j){
    cell_type_i_pca <- filter(df_PCA, cell_type == cell_type_i) |> 
      select(-cell_type)
    
    cell_type_j_pca <- filter(df_PCA, cell_type == cell_type_j) |> 
      select(-cell_type)
    
    cell_type_i_avg <- colSums(cell_type_i_pca) / nrow(cell_type_i_pca)
    cell_type_j_avg <- colSums(cell_type_j_pca) / nrow(cell_type_j_pca)
    
    wtVCosSimilarity(cell_type_i_avg, cell_type_j_avg, 
                     var_explained$`Proportion of Variance`)
  }) |> unlist()
  names(sim) <- celltypes
  sim
}) |> bind_rows()

cell_type_sims <- as.data.frame(cell_type_sims)
rownames(cell_type_sims) <- celltypes

pdf(snakemake@output$heatmap)
  Heatmap(cell_type_sims)
dev.off()

cell_type_sims |>
  rownames_to_column("cell_type1") |>
  pivot_longer(-cell_type1, values_to = "cosine_similarity", names_to = "cell_type2") |> 
  mutate(cohort = snakemake@wildcards$modality,
         seed = snakemake@wildcards$s) |>
  write_tsv(snakemake@output$tsv)


