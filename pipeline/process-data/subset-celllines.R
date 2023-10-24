library(tidyverse)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input$sce)

set.seed(1)

cell_lines <- unique(sce$CellType)[sample(1:length(unique(sce$CellType)), 10)]

sub_sce <- sce[, sce$CellType %in% cell_lines]

saveRDS(sub_sce, snakemake@output$sce)