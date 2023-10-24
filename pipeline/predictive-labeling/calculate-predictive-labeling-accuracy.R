
suppressPackageStartupMessages({
    library(tidyverse)
    library(SingleCellExperiment)
    library(yardstick)
})
source("pipeline/whatsthatcell-helpers.R")

sce <- readRDS(snakemake@input[['sce']])
if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}

if(snakemake@wildcards[['mod']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}

sce_gt <- tibble(cell_id = colnames(sce),
                 annotated_cell_type = sce$CellType)

AL_ent <- lapply(snakemake@input$AL_ent, read_tsv) |>
  bind_rows() |>
  select(-corrupted_cell_type, -iteration, -method) |>
  filter(labeling == "pred") |> 
  separate(params, c('rm_mod', 'mod', 'rm_init', 'initial', 'rm_sel', 'sel1', 'sel2', 'rm_strat', 'strat', 'rm_al', 'al',
                     'rm_pred', 'pred', 'rm_rand', 'rand', 'rm_corr', 'corrupted', 'rm_knn', 'knn', 'rm_res', 'res', 'rm_s', 'seed', 'rm_cell', 'cell_num', 'rm_pred_cell', 'pred_sel'), sep = '-') |>
  unite(selection_procedure, c(sel1, sel2), sep = "-") |>
  select(-starts_with('rm'))

AL_maxp <- lapply(snakemake@input$AL_maxp, read_tsv) |>
  bind_rows() |>
  select(-corrupted_cell_type, -iteration, -method) |>
  filter(labeling == "pred") |> 
  separate(params, c('rm_mod', 'mod', 'rm_init', 'initial', 'rm_sel', 'sel1', 'sel2', 'rm_strat', 'strat', 'rm_al', 'al',
                     'rm_pred', 'pred', 'rm_rand', 'rand', 'rm_corr', 'corrupted', 'rm_knn', 'knn', 'rm_res', 'res', 'rm_s', 'seed', 'rm_cell', 'cell_num', 'rm_pred_cell', 'pred_sel'), sep = '-') |>
  unite(selection_procedure, c(sel1, sel2), sep = "-") |>
  select(-starts_with('rm'))

random <- lapply(snakemake@input$random, read_tsv) |>
  bind_rows() |>
  select(-corrupted_cell_type, -iteration, -method) |>
  filter(labeling == "pred") |> 
  separate(params, c('rm_mod', 'mod', 'rm_init', 'initial', 'rm_sel', 'selection_procedure', 'rm_strat', 'strat', 'rm_al', 'al',
                     'rm_pred', 'pred', 'rm_rand', 'rand', 'rm_corr', 'corrupted', 'rm_knn', 'knn', 'rm_res', 'res', 'rm_s', 'seed', 'rm_cell', 'cell_num', 'rm_pred_cell', 'pred_sel'), sep = '-') |>
  select(-starts_with('rm'))

AR <- lapply(snakemake@input$AR, read_tsv) |>
  bind_rows() |>
  select(-corrupted_cell_type, -iteration, -method) |>
  filter(labeling == "pred") |> 
  separate(params, c('rm_mod', 'mod', 'rm_init', 'initial', 'rm_sel', 'sel1', 'sel2', 'rm_strat', 'strat', 'rm_al', 'al',
                     'rm_pred', 'pred', 'rm_rand', 'rand', 'rm_corr', 'corrupted', 'rm_knn', 'knn', 'rm_res', 'res', 'rm_s', 'seed', 'rm_cell', 'cell_num', 'rm_pred_cell', 'pred_sel'), sep = '-') |>
  unite(selection_procedure, c(sel1, sel2), sep = "-") |>
  select(-starts_with('rm'))
  
predictions <- bind_rows(AL_ent, AL_maxp, random, AR) |>
  dplyr::rename("predicted_cell_type" = "cell_type") |> 
  mutate(modality = snakemake@wildcards$mod) |> 
  left_join(sce_gt, by = "cell_id")

acc <- predictions %>% 
  filter(predicted_cell_type != "unassigned") %>%
  group_by(mod, initial, selection_procedure, strat, al, pred, rand, corrupted, knn, res, seed, cell_num, pred_sel) |>
  acc_wrap()

write_tsv(acc, snakemake@output$acc)