suppressPackageStartupMessages({
  library(scater)
  library(SingleCellExperiment)
  library(tidyverse)
  library(caret)
  source(file.path("utils/", 'CyTOF_LDAtrain.R'))
  source(file.path("utils/", 'CyTOF_LDApredict.R'))
})
sce_train <- readRDS(snakemake@input[['training_rds']])

labs <- read_tsv(snakemake@input[['labels']])
labs <- labs %>% column_to_rownames("cell_id")

sce_train <- sce_train[,colnames(sce_train) %in% rownames(labs)]
sce_train$cell_type <- labs[colnames(sce_train), 'cell_type']

sce_annotate <- read_rds(snakemake@input[['annotation_rds']])

data_train <- as.data.frame(t(logcounts(sce_train)))
data_train$cell_type <- sce_train$cell_type

data_annotate <- as.data.frame(t(logcounts(sce_annotate)))

train_dir <- file.path(tempdir(), "train")
annotate_dir <- file.path(tempdir(), "annotate")

dir.create(train_dir)
dir.create(annotate_dir)

write.table(data_train, file = file.path(train_dir, "train.csv"), 
            col.names = FALSE, row.names = FALSE, sep = ',')
write.table(data_annotate, file = file.path(annotate_dir, "annotate.csv"),
            col.names = FALSE, row.names = FALSE, sep = ',')

LDA.Model <- CyTOF_LDAtrain(TrainingSamplesExt = train_dir, TrainingLabelsExt = '',
                            mode = 'CSV', RelevantMarkers = seq_len(nrow(sce_train)),
                            LabelIndex = ncol(data_train), Transformation = FALSE)

predictions <- CyTOF_LDApredict(LDA.Model, TestingSamplesExt = annotate_dir,
                                mode = 'CSV', RejectionThreshold = 0)

predictions <- unlist(predictions)

if(is.null(snakemake@wildcards[['similarity']])){
  df_output <- tibble(
    cell_id = rownames(data_annotate),
    predicted_cell_type = predictions,
    prediction_params = paste0("CyTOF-LDA-knn-", snakemake@wildcards[['neighbors']], "-res-", snakemake@wildcards[['res']], 
                              '-cell_numbers-', snakemake@wildcards[['cell_num']], '-randomSelection-', snakemake@wildcards[['rand']], 
                              '-corrupted-', snakemake@wildcards[['corrupt']], '-Init-', snakemake@wildcards[['initial']], '-seed-', snakemake@wildcards[['s']]),
    selection_procedure = paste0(snakemake@wildcards[['selection_procedure']], '-strategy-', snakemake@wildcards[['strat']], '-ALAlg-', snakemake@wildcards[['AL_alg']]),
    training_annotator = snakemake@wildcards[['annotator']],
    modality = 'CyTOF'
  )
}else{
  df_output <- tibble(
    cell_id = rownames(data_annotate),
    predicted_cell_type = predictions,
    prediction_params = paste0("CyTOF-LDA-knn-", snakemake@wildcards[['neighbors']], "-res-", snakemake@wildcards[['res']], 
                              '-cell_numbers-', snakemake@wildcards[['cell_num']], '-randomSelection-', snakemake@wildcards[['rand']], 
                              '-corrupted-', snakemake@wildcards[['corrupt']], '-Init-', snakemake@wildcards[['initial']], '-seed-', snakemake@wildcards[['s']]),
    selection_procedure = paste0(snakemake@wildcards[['selection_procedure']], '-strategy-', snakemake@wildcards[['strat']], '-ALAlg-', snakemake@wildcards[['AL_alg']]),
    similarity = paste0(snakemake@wildcards[['bal']], '-', snakemake@wildcards[['similarity']]),
    training_annotator = snakemake@wildcards[['annotator']],
    modality = 'CyTOF'
  )
}

if(is.null(snakemake@wildcards[['cell_selection']])){
  df_output$cell_selection <- NA
}else{
  df_output$cell_selection <- snakemake@wildcards[['cell_selection']]
}

if(is.null(snakemake@wildcards[['similarity']]) & !is.null(snakemake@wildcards[['bal']])){
  df_output$balance <- snakemake@wildcards[['bal']]
}

if(!is.null(snakemake@wildcards$rem_percentage)){
  df_output$rem_percentage <- snakemake@wildcards$rem_percentage
}

if(!is.null(snakemake@wildcards$cell_selection)){
  df_output$pred_cells <- snakemake@wildcards$cell_selection
}

write_tsv(df_output, snakemake@output[['prediction']])