library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[['sce']])
sce_gt <- tibble(cell_id = colnames(sce),
                 cell_type = sce$CellType)

set.seed(1)
subset1 <- sce[, sample(1:ncol(sce), as.integer(snakemake@wildcards[['cell_num']]))]

subset1 <- subset1[, -c(1:2)]

each_ct <- filter(sce_gt, !(cell_id %in% colnames(subset1))) |>
  group_by(cell_type) |>
  slice_head(n = 1)

subset1 <- cbind(subset1, sce[, colnames(sce) %in% each_ct$cell_id])

gt_1 <- tibble(cell_id = colnames(subset1),
               cell_type = subset1$CellType,
               method = "random",
               set = paste0("Set", snakemake@wildcards[['set_num']]),
               params = paste0("cell_num-", snakemake@wildcards[['cell_num']], "-corruption_percentage-", snakemake@wildcards[['corrupt']]))

write_rds(subset1, snakemake@output[['sce_subset1']])
write_tsv(gt_1, snakemake@output[['gt_subset1']])