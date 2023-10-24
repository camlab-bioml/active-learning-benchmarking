suppressPackageStartupMessages({
  library(caret)
  library(SingleCellExperiment)
  library(yaml)
  library(data.table)
  library(tibble)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
})

full_pred <- read_tsv(snakemake@input$full_pred)

percent_sel <- gsub("top", "", snakemake@wildcards$cell_selection) |>
    as.numeric()
cell_num2sel <- round((percent_sel / 100) * nrow(full_pred), 0)

labelled <- filter(full_pred, labeling == "gt")

unlabelled <- filter(full_pred, labeling == "pred") |>
    arrange(entropy) |>
    slice_head(n = cell_num2sel)

bind_rows(labelled, unlabelled) |>
    mutate(corrupted_cell_type = NA, 
           iteration = NA, 
           method = NA, 
           cell_num = snakemake@wildcards$cell_num,
           params = paste0("modality-", snakemake@wildcards$modality, "-Initial-", snakemake@wildcards$initial,
                           "-selection-", snakemake@wildcards$selection_procedure, "-strat-", snakemake@wildcards$strat,
                           "-ALAlg-", snakemake@wildcards$AL_alg, "-predAlg-", snakemake@wildcards$pred_alg, 
                           "-rand_sel-", snakemake@wildcards$rand, "-corrupt-", snakemake@wildcards$corrupt,
                           "-knn-", snakemake@wildcards$neighbors, "-res-", snakemake@wildcards$res,
                           "-seed-", snakemake@wildcards$s, "-cell_num-", snakemake@wildcards$cell_num, 
                           "-cell_sel-", snakemake@wildcards$cell_selection)) |>
    write_tsv(snakemake@output$data)