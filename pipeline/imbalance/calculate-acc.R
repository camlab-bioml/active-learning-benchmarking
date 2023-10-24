suppressPackageStartupMessages({
  library(yardstick)
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

predictions <- lapply(snakemake@input$predictions, read_tsv) |>
  bind_rows() |> 
  mutate(filter_similarity = str_split_fixed(similarity, "-", 2)[,2]) |>
  filter(filter_similarity == snakemake@wildcards$similarity) |> 
  select(-filter_similarity)

sce <- readRDS(snakemake@input$sce)
if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}

if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}

cell_types <- snakemake@params |> unlist()

gt <- tibble(cell_id = colnames(sce),
             annotated_cell_type = sce$CellType) |>
  filter(annotated_cell_type %in% cell_types)

predictions_gt <- left_join(predictions, gt) |>
    separate(prediction_params, c('m1', 'm2', 'rm_knn', 'knn', 
                     'rm_res', 'res', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted', 'rm_init', 'init', 'rm_seed', 'seed'), sep = '-') |>
  unite(method, c(m1, m2), sep = '-') |> 
  mutate(selection_procedure = gsub("Active-Learning_entropy-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning_maxp-strategy-", "", selection_procedure),
         selection_procedure = gsub("-strategy-NA", "", selection_procedure),
         selection_procedure = gsub("Seurat-clustering", "Seurat_clustering", selection_procedure)) |> 
  separate(selection_procedure, c('strat', 'rm_AL', 'al'), sep = '-') |> 
  select(-starts_with('rm'))

acc <- predictions_gt |>
  filter(predicted_cell_type != "unassigned") |> 
  group_by(method, knn, res, cell_num, rand, corrupted, init, seed, strat, 
           al, cell_selection, similarity, modality) |> 
  acc_wrap()

write_tsv(acc, snakemake@output[['acc']])
