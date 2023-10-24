# Requires installation of custom fork
# devtools::install_github("Michael-Geuenich/singleCellNet")
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(readr)
  devtools::load_all("/singleCellNet")
  library(SingleCellExperiment)
})

# read in data and extract required formats
sce <- readRDS(snakemake@input$train_data)

train_labels <- read_tsv(snakemake@input[['annotation']]) %>% 
  as.data.frame() %>% 
  column_to_rownames('cell_id')

sce$cell_type <- train_labels[colnames(sce), 'cell_type']
sce <- sce[,!is.na(sce$cell_type)]
sce <- sce[,colnames(sce) %in% rownames(train_labels)]

test <- readRDS(snakemake@input$test_data)

# Extract data to singleCellNet format
SCN <- extractSCE(sce) 
stTM <- SCN$sampTab |> 
  rownames_to_column("cell_id")
expTMraw <- SCN$expDat

SCN_test <- extractSCE(test)
stTM_test <- SCN_test$sampTab |> 
  rownames_to_column("cell_id")
expTMraw_test <- SCN_test$expDat

# Train classifier
class_info <- scn_train(
  stTrain = stTM, 
  expTrain = expTMraw, 
  nTopGenes = 10, 
  nRand = 70, 
  nTrees = 1000, 
  nTopGenePairs = 25, 
  dLevel = "cell_type", 
  colName_samp = "cell_id"
)

# Predict across test set
classRes_val_all <- scn_predict(class_info[['cnProc']], expTMraw_test, nrand = 50)

# Extract predictions
classRes_val_all <- classRes_val_all |> 
  t() |> as.data.frame() |> 
  select(-rand)

# Find top predictions
cell_type_idx <- apply(classRes_val_all, 1, which.max)
names(cell_type_idx) <- NULL

cell_type_predictions <- lapply(cell_type_idx, function(x) colnames(classRes_val_all)[x]) |> 
  unlist()
classRes_val_all$cell_type <- cell_type_predictions


classRes_val_all |> 
  select(cell_type)

result <- tibble(cell_id = rownames(classRes_val_all),
                 predicted_cell_type = classRes_val_all$cell_type,
                 prediction_params = paste0("SCN-labels-knn-", snakemake@wildcards[['neighbors']],
                                            "-res-", snakemake@wildcards[['res']], 
                                            "-cell_numbers-", snakemake@wildcards[['cell_num']],
                                            '-randomSelection-', snakemake@wildcards[['rand']], 
                                            '-corrupted-', snakemake@wildcards[['corrupt']],
                                            '-Init-', snakemake@wildcards[['initial']],
                                            '-seed-', snakemake@wildcards[['s']]),
                 selection_procedure = paste0(snakemake@wildcards[['selection_procedure']],
                                              '-strategy-', snakemake@wildcards[['strat']],
                                              '-ALAlg-', snakemake@wildcards[['AL_alg']]),
                 training_annotator = snakemake@wildcards[['annotator']],
                 modality = snakemake@wildcards[['modality']])

# Needed for predictive labeling
if(is.null(snakemake@wildcards[['cell_selection']])){
  result$cell_selection <- NA
}else{
  result$cell_selection <- snakemake@wildcards[['cell_selection']]
}

if(!is.null(snakemake@wildcards[['similarity']])){
  result$similarity <- paste0(snakemake@wildcards[['bal']], '-', snakemake@wildcards[['similarity']])
}

if(is.null(snakemake@wildcards[['similarity']]) & !is.null(snakemake@wildcards[['bal']])){
  result$balance <- snakemake@wildcards[['bal']]
}

if(!is.null(snakemake@wildcards$rem_percentage)){
  result$rem_percentage <- snakemake@wildcards$rem_percentage
}

if(!is.null(snakemake@wildcards$cell_selection)){
  result$pred_cells <- snakemake@wildcards$cell_selection
}

write_tsv(result, snakemake@output[['predictions']])


