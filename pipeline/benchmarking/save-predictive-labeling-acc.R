suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(yardstick)
  library(ggplot2)
})

source("pipeline/whatsthatcell-helpers.R")

log <- file(snakemake@log[[1]], open="wt")
sink(log)

sce <- readRDS(snakemake@input[['sce']])
if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}
if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}
gt <- tibble(cell_id = colnames(sce),
             annotated_cell_type = sce$CellType)

AL_entropy <- lapply(snakemake@input[['active_learning_entropy']], read_tsv) |> 
  bind_rows() |> 
  mutate(cell_selection = "baseline")
AL_maxp <- lapply(snakemake@input[['active_learning_maxp']], read_tsv) |>
  bind_rows() |>
  mutate(cell_selection = "baseline")
AL <- bind_rows(AL_entropy, AL_maxp) |> 
  mutate(pred_lab_alg = "NA")

pred_lab <- lapply(snakemake@input[['pred_lab_files']], function(x){
  df <- read_tsv(x)
  pred_alg <- sub(".*pred_alg-", "", x)
  pred_alg <- sub("-rand_sel.*", "", pred_alg)
  df$pred_lab_alg <- pred_alg
  df
}) |> 
  bind_rows()

marker_seurat <- lapply(snakemake@input[['marker_seurat_files']], read_tsv) |> 
  bind_rows() |> 
  mutate(cell_selection = "baseline",
         pred_lab_alg = "NA")

noMarker_seurat <- lapply(snakemake@input[['NoMarker_seurat_files']], read_tsv) |> 
  bind_rows() |> 
  mutate(cell_selection = "baseline",
         pred_lab_alg = "NA")

seurat <- bind_rows(marker_seurat, noMarker_seurat)

random <- lapply(snakemake@input[['random_files']], read_tsv) |> 
  bind_rows() |> 
  mutate(cell_selection = "baseline",
         pred_lab_alg = "NA")

predictions <- bind_rows(AL, pred_lab, seurat, random) |> 
  left_join(gt) |> 
  separate(prediction_params, c('m1', 'm2', 'rm_knn', 'knn', 
                     'rm_res', 'res', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted', 'rm_init', 'initial', 'rm_s', 'seed'), sep = '-') |> 
  select(-starts_with('rm')) |> 
  unite(method, c(m1, m2), sep = '-') |> 
  mutate(selection_procedure = gsub("Active-Learning_entropy-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning_maxp-strategy-", "", selection_procedure))

acc <- predictions |> 
  filter(predicted_cell_type != "unassigned") |> 
  group_by(method, knn, res, cell_num, rand, corrupted, initial, seed, selection_procedure, 
           cell_selection, pred_lab_alg) |> 
  acc_wrap() |> 
  mutate(selection_procedure = gsub("_quant_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_quant_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("-strategy-NA", "", selection_procedure),
         cell_num = factor(cell_num, levels = c("100", "250", "500")))

write_tsv(acc, snakemake@output[['acc']])