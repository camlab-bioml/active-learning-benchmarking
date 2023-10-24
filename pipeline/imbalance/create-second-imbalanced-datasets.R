suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(as.integer(snakemake@wildcards$s))

sce <- readRDS(snakemake@input$sce)

if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}

if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}

# Cell types to select
cell_types <- snakemake@params$cell_types

# subset to cell types of interest
sce <- sce[, sce$CellType %in% cell_types]

ct_table <- tibble(cell_id = colnames(sce), cell_type = sce$CellType)

# Create balanced dataset
bal_dataset <- lapply(cell_types, function(x){
  ct <- filter(ct_table, cell_type == x)
  ct <- ct[sample(1:nrow(ct), 100), ]
  ct
}) |> bind_rows()


# Create imbalanced dataset
maj_cell_type <- snakemake@params$maj_cell_type
imbalanced_dataset <- lapply(cell_types, function(x){
  # determine how many cells to select
  if(x == maj_cell_type){ # 400 if majority
    num_cells <- 400
  }else{ # 25 if minority
    num_cells <- 25
  }
  
  # select cells
  ct <- filter(ct_table, cell_type == x)
  ct <- ct[sample(1:nrow(ct), num_cells), ]
  ct
}) |> bind_rows()


bal_sce_train <- sce[, bal_dataset$cell_id]
bal_sce_test <- sce[, !(sce$CellType %in% bal_dataset$cell_id)]
imbal_sce_train <- sce[, imbalanced_dataset$cell_id]
imbal_sce_test <- sce[, !(sce$CellType %in% imbalanced_dataset$cell_id)]

saveRDS(bal_sce_train, snakemake@output$bal_sce_train)
saveRDS(bal_sce_test, snakemake@output$bal_sce_test)
saveRDS(imbal_sce_train, snakemake@output$imbal_sce_train)
saveRDS(imbal_sce_test, snakemake@output$imbal_sce_test)

