suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(readr)
  library(dplyr)
})

set.seed(as.integer(snakemake@wildcards$s))

sce <- readRDS(snakemake@input[['sce']])

if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}

majority_ct <- snakemake@params[['majority']]
minority_ct <- snakemake@params[['minority']]

if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}

# Number of cells to select
if(snakemake@wildcards$bal == "balanced"){
  majority_num <- 250
  minority_num <- 250
}else if(snakemake@wildcards$bal == "imbalanced"){
  majority_num <- 450
  minority_num <- 50
}

# subset to majority and minority cell type datasets & sample
big_sce <- sce[, sce$CellType == majority_ct]
big_sce <- big_sce[, sample(1:ncol(big_sce), majority_num)]

small_sce <- sce[, sce$CellType == minority_ct]
small_sce <- small_sce[, sample(1:ncol(small_sce), minority_num)]

comb_sce <- cbind(big_sce, small_sce)
comb_sce <- comb_sce[, sample(1:ncol(comb_sce), ncol(comb_sce))]

# Create test set
rem_sce <- sce[, !(colnames(sce) %in% colnames(comb_sce))]
rem_sce <- rem_sce[, rem_sce$CellType == majority_ct | rem_sce$CellType == minority_ct]

table(rem_sce$CellType) |> 
  as.data.frame() |>
  dplyr::rename('cell_type' = 'Var1') |>
  mutate(params = paste0(snakemake@wildcards$modality, "-", snakemake@wildcards$bal, "-", snakemake@wildcards$similarity, "-seed-", snakemake@wildcards$s)) |>
  write_tsv(snakemake@output$test_set_tsv)

saveRDS(comb_sce, snakemake@output[['sce']])
saveRDS(rem_sce, snakemake@output[['rem_sce']])