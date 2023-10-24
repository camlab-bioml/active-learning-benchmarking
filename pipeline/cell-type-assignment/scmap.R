suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scmap)
  library(tidyverse)
  library(caret)
})

set.seed(1)

# Read in training data expression
train_sce <- readRDS(snakemake@input[['train_data']])

# Read in annotated labels
annotations <- read_tsv(snakemake@input[['annotation']]) %>% 
  as.data.frame() %>%
  column_to_rownames('cell_id')

train_sce <- train_sce[,colnames(train_sce) %in% rownames(annotations)]

# Add labels to expression
train_sce$cell_type1 <- annotations[colnames(train_sce), 'cell_type']

train_sce <- train_sce[,!is.na(train_sce$cell_type1)]

# Process data
rowData(train_sce)$feature_symbol <- rownames(train_sce)
train_sce <- selectFeatures(train_sce, suppress_plot = TRUE)

# Read in test set
annotate_sce <- readRDS(snakemake@input[['test_data']])
rowData(annotate_sce)$feature_symbol <- rownames(annotate_sce)

### Start training
train_sce <- indexCluster(train_sce)


# Map clusters
scmapCluster_results <- scmapCluster(
  projection = annotate_sce, 
  index_list = list(
    yan = metadata(train_sce)$scmap_cluster_index
  )
)

clustering_prediction <- tibble(cell_id = colnames(annotate_sce),
                                predicted_cell_type = scmapCluster_results$combined_labs,
                                prediction_params = paste0("scmap-clusters-knn-", snakemake@wildcards[['neighbors']],
                                                           "-res-", snakemake@wildcards[['res']], 
                                                           "-cell_numbers-", snakemake@wildcards[['cell_num']],
                                                           '-randomSelection-', snakemake@wildcards[['rand']], 
                                                           '-corrupted-', snakemake@wildcards[['corrupt']],
                                                           '-Init-', snakemake@wildcards[['initial']],
                                                           '-seed-', snakemake@wildcards[['s']]),
                                selection_procedure = paste0(snakemake@wildcards[['selection_procedure']],
                                                             '-strategy-', snakemake@wildcards[['strat']],
                                                             '-ALAlg-', snakemake@wildcards[['AL_alg']]),
                                modality = snakemake@wildcards[['modality']])

if(is.null(snakemake@wildcards[['cell_selection']])){
  clustering_prediction$cell_selection <- NA
}else{
  clustering_prediction$cell_selection <- snakemake@wildcards[['cell_selection']]
}

if(!is.null(snakemake@wildcards[['similarity']])){
  clustering_prediction$similarity <- paste0(snakemake@wildcards[['bal']], '-', snakemake@wildcards[['similarity']])
}

if(is.null(snakemake@wildcards[['similarity']]) & !is.null(snakemake@wildcards[['bal']])){
  clustering_prediction$balance <- snakemake@wildcards[['bal']]
}

if(!is.null(snakemake@wildcards$rem_percentage)){
  clustering_prediction$rem_percentage <- snakemake@wildcards$rem_percentage
}

if(!is.null(snakemake@wildcards$cell_selection)){
  clustering_prediction$pred_cells <- snakemake@wildcards$cell_selection
}


write_tsv(clustering_prediction, snakemake@output[['cluster_predictions']])

### Cell level
train_sce <- indexCell(train_sce)

scmapCell_results <- scmapCell(
  annotate_sce, 
  list(
    yan = metadata(train_sce)$scmap_cell_index
  )
)

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(train_sce)$cell_type1)
  )
)

sc_prediction <- tibble(cell_id = colnames(annotate_sce),
                        predicted_cell_type = scmapCell_clusters$scmap_cluster_labs[,1],
                        prediction_params = paste0("scmap-sc-knn-", snakemake@wildcards[['neighbors']],
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

if(is.null(snakemake@wildcards[['cell_selection']])){
  sc_prediction$cell_selection <- NA
}else{
  sc_prediction$cell_selection <- snakemake@wildcards[['cell_selection']]
}

if(!is.null(snakemake@wildcards[['similarity']])){
  sc_prediction$similarity <- paste0(snakemake@wildcards[['bal']], '-', snakemake@wildcards[['similarity']])
}

if(is.null(snakemake@wildcards[['similarity']]) & !is.null(snakemake@wildcards[['bal']])){
  sc_prediction$balance <- snakemake@wildcards[['bal']]
}

if(!is.null(snakemake@wildcards$rem_percentage)){
  sc_prediction$rem_percentage <- snakemake@wildcards$rem_percentage
}

if(!is.null(snakemake@wildcards$cell_selection)){
  sc_prediction$pred_cells <- snakemake@wildcards$cell_selection
}

write_tsv(sc_prediction, snakemake@output[['sc_predictions']])

