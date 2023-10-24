suppressPackageStartupMessages({
  library(tidyverse)
  library(yardstick)
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

sce_gt <- tibble(cell_id = colnames(sce),
                 annotated_cell_type = sce$CellType)

predictions <- lapply(snakemake@input[['predictions']], read_tsv) |>
  bind_rows()

predictions <- predictions |> 
  left_join(sce_gt) %>% 
  separate(prediction_params, c('m1', 'm2', 'rm_knn', 'knn', 
                     'rm_res', 'res', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted', 'rm_init', 'initial', 'rm_seed', 'seed'), sep = '-') %>% 
  select(-starts_with('rm')) %>% 
  unite(method, c(m1, m2), sep = '-') %>% 
  mutate(selection_procedure = gsub("Active-Learning_entropy-strategy-", "", selection_procedure),
         selection_procedure = gsub("Active-Learning_maxp-strategy-", "", selection_procedure),
         selection_procedure = gsub("Seurat-clustering", "Seurat_clustering", selection_procedure), #|> 
         selection_procedure = gsub("-strategy-NA", "", selection_procedure)) |>
  separate(selection_procedure, c('selection_procedure', 'rm_alg', 'AL_alg'), sep = '-') |>
  select(-starts_with('rm'))

acc <- predictions %>% 
  filter(predicted_cell_type != "unassigned") %>%
  group_by(method, knn, res, cell_num, rand, corrupted, initial, seed, selection_procedure, AL_alg) %>% 
  acc_wrap() %>% 
  mutate(selection_procedure = gsub("_quant_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_entropy", "-entropy-AL", selection_procedure),
         selection_procedure = gsub("_quant_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("_maxp", "-maxp-AL", selection_procedure),
         selection_procedure = gsub("-strategy-NA", "", selection_procedure),
         cell_num = factor(cell_num, levels = c("100", "250", "500")))

write_tsv(acc, snakemake@output[['acc']])