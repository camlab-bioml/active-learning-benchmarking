suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(caret)
})

sce <- readRDS(snakemake@input[['rds']])

if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}

if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}

set.seed(as.integer(snakemake@wildcards[['s']]))
sce <- sce[, sample(1:ncol(sce))]
train <- createDataPartition(sce$CellType, p = 0.5)$Resample1

train_sce <- sce[,train]
test_sce <- sce[,-train]

saveRDS(train_sce, snakemake@output[['train']])
saveRDS(test_sce, snakemake@output[['test']])