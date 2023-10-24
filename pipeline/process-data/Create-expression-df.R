library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[['expression_sce']])

expression <- assays(sce)$logcounts %>% 
  t() %>% as.matrix() %>% 
  as.data.frame()

expression$cell_type <- sce[, rownames(expression)]$CellType
expression$cell_id <- rownames(expression)
rownames(expression) <- NULL

write_tsv(expression, snakemake@output[['expression_df']])